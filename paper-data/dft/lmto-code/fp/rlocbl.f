C#define F90
      subroutine rlocbl(s_site,s_spec,s_lat,lfrce,nbas,isp,q,ndham,
     .  ndimh,nspc,napw,igvapw,iprmb,numq,nevec,evec,ewgt,evl,lcplxp,
     .  lekkl,f)
C- Accumulates the local atomic densities.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos qkkl qhkl qhhl eqkkl eqhkl eqhhl pikk sigkk
Ci                 pihk sighk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa kmxt lmxb rsma pz name orbp
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
Ci Inputs
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   nbas  :size of basis
Ci   isp   :spin channel
Ci   q     :Bloch wave number
Ci   ndham :leading dimension of evl
Ci   ndimh :dimension of evec
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial fermi levels
Ci   nevec :number of occupied eigenvectors
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   evl   :eigenvalues
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lekkl :0 do not accumulate eqkkl; 1 do accumulate eqkkl
Co Outputs
Co   f     :local contribution to forces is added
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :spin index in both independent and coupled spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cl         :   Collinear case        Noncollinear case
Cl         :   isp = spin index      isp  = 1
Cl         :   ispc = 1              ispc = spin index
Cl         :   ksp = spin index      ksp  = spin index
Cr Remarks
Cr   The density is constructed from eigenstates as sum_n^occ [psi_n psi+_n]
Cr   In the LMTO basis it is represented through the density matrix
Cr      Dij = {sum_n w_n evec_in [evec_jn]+}
Cr   This routine makes qkkl:  they are contractions of Dij and
Cr   the coefficients to the one-center expansion of the wave
Cr   function inside the augmentation sphere:
Cr     F~i = Fi + sum_kL C^i_kL (P~kL - PkL)
Cr   As usual, we neglect cross terms when making function products.
Cr   Thus function products are of the form
Cr     F~i F~j = Fi Fj +
Cr             = sum_kLk'L' C^i_kL (P~kL P~k'L' - PkL Pk'L') C^j_k'L'
Cr             = sum_kLk'L' C^i_kL (n1kLk'L' - n2kLk'L') C^j_k'L'
Cr   the qkkl are defined as, e.g.
Cr      qpp_kLk'L' = sum_ij D_ij C^i_kL C+^j_k'L'
Cr   so that the local part of the output density is
Cr      n1 - n2 = sum_kLk'L' qpp_kLk'L' (n1kLk'L' - n2kLk'L')
Cu Updates
Cu   17 Dec 13 Write qkkl as C^i_kL C+^j_k'L' instead of C+^i_kL C^j_k'L'
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   22 Nov 12 Replace qkkl with structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   09 Apr 10 First attempt at noncollinear smrho.
Cu             Not checked; force terms not implemented
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   16 Jun 05 Makes spin-off-diagonal density matrix, noncollinear case
Cu   23 Dec 04 Extended to spin-coupled case
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 Jun 00 spin polarized
Cu   25 May 00 Adapted from nfp rloc_q.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPI
C      include "mpif.h"
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      integer procid,master
C      integer numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1), pid
C      double precision starttime, endtime
C      character*120 strn
C      logical mlog,cmdopt
C#endif
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
      integer lfrce,nbas,isp,ndimh,nspc,numq,nevec,lcplxp,lekkl,
     .  iprmb(1),ndham,napw,igvapw(3,napw)
C     integer osig(3,1),otau(3,1),oppi(3,1)
      double precision ewgt(numq,nevec),evl(ndham,isp),f(3,nbas,numq),
     .  q(3)
      double complex evec(ndimh,nspc,nevec)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(nbas)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat

C ... Dynamically allocated local arrays
      real(8), pointer :: qkk(:,:),qhk(:,:),qhh(:,:)
      real(8), pointer :: eqkk(:,:),eqhk(:,:),eqhh(:,:)
      real(8), pointer :: ppikk(:,:),ppihk(:,:) !,ppihh(:,:)
      real(8), pointer :: sigkk(:,:),sighk(:,:) !,sighh(:,:)
      complex(8), allocatable :: b(:)
      complex(8), allocatable :: db(:)
      complex(8),allocatable:: cPkL(:),da(:),wk(:,:)

C ... Local parameters or process-shared variables
C#ifdefC MPI
C      integer, dimension(:), allocatable :: bproc
C      integer nelt(3),lgunit,lmxh,nlmh,nsp,i,k,is
C#endif
      integer nlmbx,nlmx,ktop0,npmx,nkap0,n0
      parameter (nlmbx=25, npmx=32, nkap0=4, n0=10)
      integer kmaxx,nlmax,nglob,nlmto
      double precision alat,qlat(3,3),xx
      double precision, allocatable :: floc(:)
C ... Local process-specific variables
      integer ia,isa,ivec,kmax,lmxa,nlma,
     .        lmxha,nlmha,nkaph,ispc,ksp
      double precision pa(3),rsma,pi,tpiba

C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,strn)
C      nsp = nglob('nsp')
C      if (mlog) then
C        do  pid = 0, numprocs-1
C          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,
C     .                   MPI_COMM_WORLD,ierr)
C          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
C     .                   MPI_COMM_WORLD,ierr)
C        enddo
C      endif
C#endif

      if (nevec <= 0) return
      call tcn('rlocbl')

C --- Setup ---
      nkaph = nglob('nkaph')

C ... Find maximum sizes needed to allocate strux; allocate them
      call lkdim(1,nbas,s_site,s_spec,kmaxx,nlmax)
      nlmax = (nlmax+1)**2
C      nlmax = 0
C      kmaxx = 0
C      do  ia = 1, nbas
C        isa = s_site(ia)%spec
C        lmxa = s_spec(isa)%lmxa
C        kmax = s_spec(isa)%kmxt
C        nlma = (lmxa+1)**2
C        kmaxx = max(kmaxx,kmax)
C        nlmax = max(nlmax,nlma)
C      enddo
      nlmto = ndimh-napw
C     Needed for PW part
      alat = s_lat%alat
      qlat = s_lat%qlat
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat

C ... Allocate workspace for augmentation arrays
      nlmx  = nlmax
      ktop0 = kmaxx
C      allocate(cPkL(0:ktop0,nlmx,nspc),da(0:ktop0,nlmx,3),
      allocate(cPkL((ktop0+1)*nlmx*nspc),da((ktop0+1)*nlmx*3),
     .  wk(0:ktop0,nlmx))
      if (nlmax > nlmx) call rxi('rlocbl: nlmx < nlma=',nlmax)
      if (kmaxx > ktop0) call rxi('rlocbl: ktop0 < kmax=',kmax)

C ... Allocate workspace for strux
      allocate(b(ndimh*nlmax*(kmaxx+1)))
      allocate(db(ndimh*nlmax*(kmaxx+1)*3))
      if (lfrce /= 0) then
        allocate(floc(3*nbas*numq))
        call dpzero(floc,3*nbas*numq)
      endif

C --- Loop over augmentation sites ---
C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_RLOCBL,procid,"rlocbl")
C#endifC
C      allocate (bproc(0:numprocs), stat=ierr)
C      call dstrbp(nbas,numprocs,1,bproc(0))
C      do  ia = bproc(procid), bproc(procid+1)-1
C        if (mlog .and. ia == bproc(procid)) then
C          call gettime(datim)
C          call awrit4(' rlocbl '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' starting atoms %i to %i',' ',256,lgunit(3),
C     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
C        endif
C#else
      do  ia = 1, nbas
C#endif

        isa = s_site(ia)%spec
        pa = s_site(ia)%pos
        lmxa = s_spec(isa)%lmxa
        if (lmxa == -1) cycle
        lmxha = s_spec(isa)%lmxb
        kmax = s_spec(isa)%kmxt
        rsma = s_spec(isa)%rsma
        nlmha = (lmxha+1)**2
        nlma  = (lmxa+1)**2

        qkk => s_site(ia)%qkkl
        qhk => s_site(ia)%qhkl
        qhh => s_site(ia)%qhhl
        eqkk => s_site(ia)%eqkkl
        eqhk => s_site(ia)%eqhkl
        eqhh => s_site(ia)%eqhhl

        ppikk => s_site(ia)%pikk
        sigkk => s_site(ia)%sigkk
        ppihk => s_site(ia)%pihk
        sighk => s_site(ia)%sighk

C   --- Strux to expand all orbitals and their gradients at site ia ---
        call bstrux(1,s_lat,s_site,s_spec,s_lat%cg,s_lat%indxcg,
     .    s_lat%jcg,s_lat%cy,iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,
     .    napw,igvapw,b,db)
C       if (ia == 2) call snott

C   --- Loop over eigenstates ---
        do  ivec = 1, nevec

C     ... Pkl expansion of eigenvector
          call rlocb1(ndimh,nlma,kmax,evec(1,1,ivec),nspc,b,cPkL)

C         In noncollinear case, isp=1 always => need internal ispc=1..2
C         See local variable definitions above for ksp
          do  ispc = 1, nspc
          ksp = max(ispc,isp)

          if (ispc == 2)
     .    call dswap(2*(1+kmax)*nlmx,cPkL(1),1,cPkL(1+(kmax+1)*nlma),1)

C     ... Add to local density coefficients for one state
          call prlcb3(0,kmax,nlma,ksp,nspc,cPkL,numq,ewgt(1,ivec),xx,qkk)
          call prlcb2(0,ia,nkaph,iprmb,nlmha,kmax,nlma,ksp,nspc,cPkL,
     .      nlmto,ndimh,evec(1,ispc,ivec),ewgt(1,ivec),numq,xx,qhh,qhk)
          if (lekkl == 1) then
            call prlcb3(1,kmax,nlma,ksp,1,cPkL,numq,ewgt(1,ivec),
     .        evl(ivec,isp),eqkk)
            call prlcb2(1,ia,nkaph,iprmb,nlmha,kmax,nlma,ksp,1,cPkL,
     .        nlmto,ndimh,evec(1,ispc,ivec),ewgt(1,ivec),numq,
     .        evl(ivec,isp),eqhh,eqhk)
          endif

C     ... Contribution to forces
          if (lfrce /= 0) then
            call rxx(nspc /= 1,'forces not implemented in noncoll case')
            call flocbl(nbas,ia,kmax,nkaph,lmxha,nlmha,nlma,lmxa,nlmto,
     .        ndimh,iprmb,ksp,evl(ivec,isp),evec(1,ispc,ivec),
     .        ewgt(1,ivec),numq,cPkL,db,da,wk,ppikk,ppikk,
     .        sigkk,ppihk,ppihk,sighk,lcplxp,floc)
          endif
C         print *, 'ia',ia,ivec,ispc,sngl(qkk(1,1:2))
        enddo
        enddo
      enddo  ! loop over ia
      deallocate(b,db)
C      ia = size(qhk)/4
C      call prmx('qhk(1)',qhk(1,1),ia,ia,1)
C      call prmx('qhk(2)',qhk(1,2),ia,ia,1)
C      call prmx('qhk(3)',qhk(1,3),ia,ia,1)
C      call prmx('qhk(4)',qhk(1,4),ia,ia,1)
C      call prmx('qhk(5)',qhk(1,5),ia,ia,1)
C      call prmx('qhk(6)',qhk(1,6),ia,ia,1)
C      stop

C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_RLOCBL,procid,"rlocbl")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endifC
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_BCAST,procid,"broadcast")
C#endifC
C      do  pid = 0, numprocs-1
C        do  ia = bproc(pid), bproc(pid+1)-1
C          is = s_site(ia)%spec
C          lmxa = s_spec(is)%lmxa
C          if (lmxa == -1) cycle
C          lmxh = s_spec(is)%lmxb
C          kmax = s_spec(is)%kmxt
C          nlma = (lmxa+1)**2
C          nlmh = (lmxh+1)**2
C          nelt(1) = (kmax+1)*(kmax+1)*nlma*nlma
C          nelt(2) = (kmax+1)*nkaph*nlma*nlmh
C          nelt(3) = nkaph*nkaph*nlmh*nlmh
C          i = numq*nsp*nspc
C          k = numq*nsp*(2*nspc-1)
C          call mpibc3(s_site(ia)%qkkl,nelt(1)*i,4,pid,.false.,' ',' ')
C          call mpibc3(s_site(ia)%qhkl,nelt(2)*k,4,pid,.false.,' ',' ')
C          call mpibc3(s_site(ia)%qhhl,nelt(3)*i,4,pid,.false.,' ',' ')
C          if (lekkl == 1) then
C           call mpibc3(s_site(ia)%eqkkl,nelt(1)*i,4,pid,.false.,' ',' ')
C           call mpibc3(s_site(ia)%eqhkl,nelt(2)*i,4,pid,.false.,' ',' ')
C           call mpibc3(s_site(ia)%eqhhl,nelt(3)*i,4,pid,.false.,' ',' ')
C          endif
C        enddo
C      enddo
C      if (lfrce /= 0) then
C        call mpibc2(floc,3*nbas*numq,4,3,.false.,' ',' ')
CC        allocate(buffer(3*nbas*numq), stat=ierr)
CC        call MPI_ALLREDUCE(floc,buffer,3*nbas*numq,
CC     .       mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
CC        if (mlog) then
CC          call gettime(datim)
CC          call awrit2(' rlocbl '//datim//' Process %i of %i on '
CC     .        //shortname(procid)(1:namelen(procid))//
CC     .        ' allreduce forces',' ',256,lgunit(3),
CC     .        procid,numprocs)
CC        endif
C        call daxpy(3*nbas*numq,1d0,floc,1,f,1)
C        deallocate(floc, stat=ierr)
C      endif
C      deallocate(bproc, stat=ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BCAST,procid,"broadcast")
C#endifC
C#else
      if (lfrce /= 0) then
        call daxpy(3*nbas*numq,1d0,floc,1,f,1)
          deallocate(floc)
      endif
C#endif

      deallocate(cPkL,da,wk)
      call tcx('rlocbl')

      end

      subroutine rlocb1(ndimh,nlma,kmax,evec,nspc,b,cPkL)
C- Pkl expansion of wave function at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimh :dimension of evec
Ci   nlma  :augmentation L-cutoff in PkL expansion
Ci   kmax  :k- cutoff in PkL expansion
Ci   evec  :eigenvector coefficients
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   b     :strux to expand of orbitals from other sites in PkL
Ci         :b = b(ndimh,nlma,0:kmax)
Co Outputs
Co   cPkL  :coefficients to PkL expansion of eigenvector
Cu Updates
Cu   09 Apr 10 Returns cPkL for both spins in noncollinear case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,ndimh,nlma,nspc
      double complex b(ndimh,nlma,0:kmax),cPkL(0:kmax,nlma,nspc),
     .  evec(ndimh,nspc)
C ... Local parameters
      integer i,k,ilma,ispc


      call tcn('rlocb1')
C     call zprm('b',2,b,ndimh,ndimh,nlma*(kmax+1))
      call dpzero(cPkL, 2*(kmax+1)*nlma*nspc)
      do  ispc = 1, nspc
      do  k = 0, kmax
        do  ilma = 1, nlma
          do  i = 1, ndimh
            cPkL(k,ilma,ispc) = cPkL(k,ilma,ispc) +
     .                          evec(i,ispc)*b(i,ilma,k)
          enddo
        enddo
      enddo
      enddo
C     call zprm('cPkL',2,cPkL,kmax+1,kmax+1,nlma*nspc)
      call tcx('rlocb1')
      end

      subroutine prlcb2(job,ia,nkaph,iprmb,nlmha,kmax,nlma,ksp,nspc,
     .  cPkL,nlmto,ndimh,evec,ewgt,numq,evl,qhh,qhp)
C- Add one and two-center terms to density coeffs
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   ia    :site of augmentation
Ci   nkaph :dimensions qhh,qhp
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions qhh,qhp
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   ksp   :spin index in both collinear and noncollinear cases
Ci         :   Collinear case        Noncollinear case
Ci         :   isp = spin index      isp  = 1
Ci         :   ispc = 1              ispc = spin index
Ci         :   ksp = spin index      ksp  = spin index
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   cPkL  :PkL expansion eigenvector at site ia.
Ci   nlmto :dimension of lmto component of basis
Ci   evec  :eigenvector
Ci   numq  :number of trial fermi levels
Ci   ewgt  :eigenvector weight
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qhh   :one-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Co         :In the noncollinear case,
Co         : qhh(:,:,:,:,:,3) = Re(qhh) spin block (1,2)
Co         : qhh(:,:,:,:,:,4) = Im(qhh) spin block (1,2)
Co         : Note: qhh(2,1) = qhh(1,2)*
Co   qhp   :two-center density-matrix for PkL expansion (job=0)
Co         :energy-weighted matrix (job=1)
Co         :qhp = qhp(ik=1:nkap,k=0:kmax,ilm=1:nlmb,ilma=1:nlma,iq,isp)
Co         :Collinear case: for each spin and state n with occ weight wn,
Co         :qhp(ik,k,ilm,ilma) += z(ik,ilm) C+(k,ilma) wn
Co         :Noncollinear case:
Co         :qhp(ik,k,ilm,ilma;1,2) += z(ik,ilm,1) C+(k,ilma,2) wn
Co         :qhp(ik,k,ilm,ilma;2,1) += z(ik,ilm,2) C+(k,ilma,1) wn
Co         :Store as:
Co         : qhp(:,:,:,:,:,3) =  Re(qhp(1,2))
Co         : qhp(:,:,:,:,:,4) =  Im(qhp(1,2))
Co         : qhp(:,:,:,:,:,5) =  Re(qhp(2,1))
Co         : qhp(:,:,:,:,:,6) = -Im(qhp(2,1))
Cr Remarks
Cr   In the collinear case, loop  ksp = 1, nsp
Cr   In the nocollinear case, loop ksp = 1, nspc
Cu Updates
Cu   17 Dec 13 Extended to noncollinear case
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
Cu   27 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,ia,kmax,nkaph,ksp,nspc,nlmto,ndimh,nlma,nlmha,numq,
     .  iprmb(1)
      double precision qhh(nkaph,nkaph,nlmha,nlmha,numq,4),
     .  qhp(nkaph,0:kmax,nlmha,nlma,numq,6),ewgt(numq),evl
      double complex evec(ndimh,nspc),cPkL(0:kmax,nlma,nspc)
C ... Local parameters
      integer i1,i2,ilm1,ilm2,ilma,io1,io2,iq,k,ik1,ik2,
     .  l1,l2,n0,nkap0,nlm11,nlm12,nlm21,nlm22
      parameter (n0=10,nkap0=4)
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .  blks(n0*nkap0),ntab(n0*nkap0)
      double precision xx
      double complex cxx

      if (nlmto == 0) return
      call tcn('prlcb2')

C --- Loop over all orbitals centered at site ia, incl. local orbs ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

      do  io1 = 1, norb

        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
C       i1 = hamiltonian offset for first orbital in block
        i1 = offl(io1)
        do  ilm1 = nlm11, nlm12
          i1 = i1+1
C     ... Accumulate products H*Pkl
C         Origin of extra factor of 2: w.f. expressed as sum (head + tail)
C         Products contain (head*head + tail*tail + 2*head*tail)
          do  iq = 1, numq
            if (job == 0) then
              do  k = 0, kmax
                do  ilma = 1, nlma
                  qhp(ik1,k,ilm1,ilma,iq,ksp) =
     .            qhp(ik1,k,ilm1,ilma,iq,ksp) +
     .              2*evec(i1,1)*dconjg(cPkL(k,ilma,1))*ewgt(iq)
                enddo
              enddo
            else
              do  k = 0, kmax
                qhp(ik1,k,ilm1,ilm1,iq,ksp) =
     .          qhp(ik1,k,ilm1,ilm1,iq,ksp) +
     .            2*evl*evec(i1,1)*dconjg(cPkL(k,ilm1,1))*ewgt(iq)
              enddo
            endif
C           Spin 12 and spin 21 blocks
            if (nspc == 2 .and. job == 0 .and. ksp == 1) then
              do  k = 0, kmax
                do  ilma = 1, nlma
                  cxx = 2*evec(i1,1)*dconjg(cPkL(k,ilma,2))*ewgt(iq)
                  qhp(ik1,k,ilm1,ilma,iq,3) =
     .            qhp(ik1,k,ilm1,ilma,iq,3) + dble(cxx)
                  qhp(ik1,k,ilm1,ilma,iq,4) =
     .            qhp(ik1,k,ilm1,ilma,iq,4) + dimag(cxx)
                  cxx = 2*evec(i1,2)*dconjg(cPkL(k,ilma,1))*ewgt(iq)
                  qhp(ik1,k,ilm1,ilma,iq,5) =
     .            qhp(ik1,k,ilm1,ilma,iq,5) + dble(cxx)
                  qhp(ik1,k,ilm1,ilma,iq,6) =
     .            qhp(ik1,k,ilm1,ilma,iq,6) - dimag(cxx)
                enddo
              enddo
            endif
          enddo

C     ... Accumulate products H*H
          do  io2 = 1, norb

            l2  = ltab(io2)
            ik2 = ktab(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks(io2)-1
C           i2 = orbital index in iprmb order
            i2 = offl(io2)
            do  ilm2 = nlm21, nlm22
              i2 = i2+1
              if (job == 0) then
                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,iq,ksp) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,ksp) +
     .              evec(i1,1)*dconjg(evec(i2,1))*ewgt(iq)
                enddo
              elseif (job == 1 .and. ilm1 == ilm2) then
                do  iq = 1, numq
                  qhh(ik1,ik2,ilm1,ilm2,iq,ksp) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,ksp) +
     .              evl*evec(i1,1)*dconjg(evec(i2,1))*ewgt(iq)
                enddo
              endif
C             Spin off-diagonal part
              if (nspc == 2 .and. job == 0 .and. ksp == 1) then
                do  iq = 1, numq
                  cxx = evec(i1,1)*dconjg(evec(i2,2))*ewgt(iq)
                  qhh(ik1,ik2,ilm1,ilm2,iq,3) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,3) + dble(cxx)
                  qhh(ik1,ik2,ilm1,ilm2,iq,4) =
     .            qhh(ik1,ik2,ilm1,ilm2,iq,4) + dimag(cxx)
                enddo
              endif
            enddo
          enddo

        enddo
      enddo

      call tcx('prlcb2')
      end

      subroutine prlcb3(job,kmax,nlma,ksp,nspc,cPkL,numq,ewgt,evl,qpp)
C- Add to local density coefficients for one state
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0 accumulate local density-matrix
Ci         :1 accumulate local density-matrix weighted by energy
Ci   kmax  :polynomial cutoff in PkL expansion
Ci   nlma  :L cutoff in PkL expansion
Ci   ksp   :spin index in both collinear and noncollinear cases
Ci         :   Collinear case        Noncollinear case
Ci         :   isp = spin index      isp  = 1
Ci         :   ispc = 1              ispc = spin index
Ci         :   ksp = spin index      ksp  = spin index
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   cPkL  :coefficients to PkL expansion of evec
Ci   numq  :number of trial fermi levels
Ci   ewgt  :eigenvector weights
Ci   evl   :energy weight (job=1)
Co Outputs
Co   qpp   :local density matrix for PkL expansion (job=0)
Co         :energy-weighted local density matrix (job=1)
Co         :In the noncollinear case,
Co         : qpp(:,:,:,:,:,3) = Re(qpp(1,2)) spin block
Co         : qpp(:,:,:,:,:,4) = Im(qpp(1,2)) spin block
Co         : Note: qpp(2,1) = qpp(1,2)*
Cr Remarks
Cu Updates
Cu   17 Dec 13 Extended to noncollinear case
Cu   05 Jul 08 (T. Kotani)
Cu             Option to accumulate energy-weighted output density
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer job,kmax,nlma,ksp,nspc,numq
      double complex cPkL(0:kmax,nlma,nspc)
      double precision qpp(0:kmax,0:kmax,nlma,nlma,numq,4),
     .  ewgt(numq),evl
C ... Local parameters
      double complex cxx
      double precision fac
      integer iq,ilm2,ilm1,k1,k2

      call tcn('prlcb3')
      do  iq = 1, numq
        fac = ewgt(iq)
        if (job == 1) fac = evl*ewgt(iq)
        do  ilm2 = 1, nlma
          do  ilm1 = 1, nlma
            do  k1 = 0, kmax
              do  k2 = 0, kmax
                qpp(k1,k2,ilm1,ilm2,iq,ksp)= qpp(k1,k2,ilm1,ilm2,iq,ksp)
     .            + fac*cPkL(k1,ilm1,1)*dconjg(cPkL(k2,ilm2,1))
              enddo
            enddo
          enddo
        enddo
        if (nspc == 2 .and. job == 0 .and. ksp == 1) then
          do  ilm2 = 1, nlma
          do  ilm1 = 1, nlma
            do  k1 = 0, kmax
            do  k2 = 0, kmax
              cxx = fac*cPkL(k1,ilm1,ksp)*dconjg(cPkL(k2,ilm2,3-ksp))
              qpp(k1,k2,ilm1,ilm2,iq,ksp+2) =
     .        qpp(k1,k2,ilm1,ilm2,iq,ksp+2) + dble(cxx)
              qpp(k1,k2,ilm1,ilm2,iq,ksp+3) =
     .        qpp(k1,k2,ilm1,ilm2,iq,ksp+3) + dimag(cxx)
              enddo
            enddo
          enddo
          enddo
        endif
      enddo
      call tcx('prlcb3')

      end

      subroutine lkdim(ib1,ib2,s_site,s_spec,kmxax,lmxax)
C- Return largest kmax and lxma in a range of sites.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa kmxt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   ib1   :Return maxval of lmxa,kmxa for sites ib1:ib2
Ci   ib2   :Return maxval of lmxa,kmxa for sites ib1:ib2
Co Outputs
Co   kmxax :maxval of kmxa
Co   lmxax :maxval of lmxa
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   31 Dec 13
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ib1,ib2,kmxax,lmxax
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(ib2)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,lmxa,kmxa

      kmxax = -1; lmxax = -1
      do  ib = ib1, ib2
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        if (lmxa == -1) cycle
        kmxa = s_spec(is)%kmxt
        kmxax = max(kmxax,kmxa)
        lmxax = max(lmxa,lmxax)
      enddo

      end
