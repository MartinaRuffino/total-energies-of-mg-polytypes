      subroutine fsmbpw(nbas,s_site,s_spec,vavg,ndimh,nlmto,iprmb,numq,
     .  nevec,evl,evec,ewgt,napw,qpgv,qpg2v,ylv,nlmax,lmxax,alat,sqv,f)
C- Force from smoothed hamiltonian (constant potential), PW contribution
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Ci Inputs
Ci   nbas  :size of basis
Ci   vavg  :constant potential (MT zero) to be added to h
Ci   ndimh :dimension of hamiltonian
Ci   nlmto :dimension of lmto part of hamiltonian
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial Fermi levels
Ci   nevec :number of occupied eigenvectors
Ci   evl   :eigenvalues
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   napw  :number of augmented PWs in basis
Ci   qpgv  :q+G
Ci   qpg2v
Ci   ylv
Ci   nlmax
Ci   lmxax
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   sqv   :square root of volume
Co Outputs
Co   f     :PW contribution to force is added to f
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   04 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPI
C      include "mpif.h"
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      integer pid, procid, master, numprocs, ierr,
C     .        status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen,i,lgunit
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C      double precision starttime, endtime
C      logical mlog,cmdopt
C      character*120 strn
C#endif
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
      integer nbas,ndimh,napw,nlmax,nlmto,nevec,numq,iprmb(nlmto),lmxax
      double precision ylv(napw,nlmax)
      double precision evl(ndimh),f(3,nbas,numq),ewgt(numq,nevec),vavg
      double complex evec(ndimh,ndimh)
      double precision qpgv(3,napw),qpg2v(napw),qpg2,alat,sqv
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer n0,nkap0
      parameter (n0=10, nkap0=4)
      integer i1,i2,ib1,ilm1,io1,iq,is1,l1,ik1,ig,ivec,nglob,
     .  nlm11,nlm12,m
      integer norb1,ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),
     .  blks1(n0*nkap0),ntab1(n0*nkap0),lh1(nkap0),nkap1
      double precision gam,denom,pi,fpi,ddot
      double precision e1(n0,nkap0),rsm1(n0,nkap0),p1(3),xx(n0),wt
      double complex phase,srm1,fach,ovl,ccc(3),sum,srm1l(0:n0)
      parameter (srm1=(0d0,1d0))
C#ifdefC MPI
C      double precision, dimension(:,:,:), allocatable :: xf
C      double precision, dimension(:),     allocatable :: buffer
C      integer, dimension(:,:), allocatable :: index
C      integer iloop,ib
C#endif

C --- Setup ---
      if (nevec <= 0) return
      call tcn ('fsmbpw')

C#ifdefC MPI
C      allocate(xf(1:3,1:nbas,1:numq), stat=ierr)
C      call dcopy(3*nbas*numq,0d0,0,xf,1)
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,strn)
C      if (mlog) then
C        do  pid = 0, numprocs-1
C          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,
C     .                   MPI_COMM_WORLD,ierr)
C          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
C     .                   MPI_COMM_WORLD,ierr)
C        enddo
C      endif
C#endif

      nbas  = nglob('nbas')
      pi = 4d0*datan(1d0)
C     tpiba = 2d0*pi/alat
      fpi = 4*pi
      srm1l(0) = 1d0
      do  l1 = 1, lmxax
        srm1l(l1) = (srm1)**l1
      enddo

C --- Loop over first and second site indices ---
C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_FSMBPW,procid,"fsmbpw")
C#endifC
C      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
C      call dstrbp(nbas,numprocs,-1,index(0,0))
C      do  iloop = 1, index(procid,0)
C        ib1 = index(procid,iloop)
C        if (mlog) then
C          call gettime(datim)
C          call awrit4(' fsmbpw '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' starting atom %i of %i',' ',256,lgunit(3),
C     .        procid,numprocs,ib1,index(procid,0))
C        endif
C#else
      do  ib1 = 1, nbas
C#endif
        is1 = s_site(ib1)%spec
        p1 = s_site(ib1)%pos
        call uspecb(0,2,s_spec,is1,is1,lh1,rsm1,e1,nkap1)
C       Row info telling fsmbpw where to poke s0 made by hhibl
        call orbl(ib1,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
C       For now, do not allow l blocks to be grouped.
C       To do so will require rewriting loops below.
        call gtbsl1(8+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)

C   ... Hsm (i1) \times i(q+G)[(q+G)**2+const] PW (i2) Takao. Taken from smhsbl.f
C       i1--> Hsm, i2--> PW
        do  ig = 1, napw
          i2 = ig + nlmto
          qpg2 = qpg2v(ig)
          phase = exp(srm1*alat*ddot(3,qpgv(1,ig),1,p1,1))
          do  io1 = 1, norb1
          if (blks1(io1) /= 0) then
C           l1,ik1 = l and kaph indices, needed to locate block in s0
            l1  = ltab1(io1)
            ik1 = ktab1(io1)
            nlm11 = l1**2+1
            nlm12 = nlm11 + blks1(io1)-1
C           Note:  using srm1l => l must be fixed in ilm loop below
            denom = e1(l1+1,ik1) - qpg2
            gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
            fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
            i1 = offl1(io1)
            do  ilm1 = nlm11, nlm12
              i1 = i1+1
C             s(i1,i2) = ovl
C             h(i1,i2) = (qpg2 + vavg) * ovl
              ovl = fach * ylv(ig,ilm1)/sqv ! Eq. 9.4 in JMP39 3393

C         ... Loop over occupied eigenstates and x,y,z
              do  m = 1, 3
              sum = 0
              do  ivec = 1, nevec

C             gradient PW * (H - E S)
              ccc(m) = ovl * srm1*qpgv(m,ig) * (qpg2 + vavg - evl(ivec))
              sum = dconjg(evec(i1,ivec))*ccc(m)*evec(i2,ivec)

              do  iq = 1, numq
                wt = ewgt(iq,ivec)
C#ifdefC MPI
C                      xf(m,ib1,iq) = xf(m,ib1,iq) - 2*wt*sum
C#else
                      f(m,ib1,iq) = f(m,ib1,iq) - 2*wt*sum
C#endif
              enddo

            enddo               !ivec
            enddo               !m
          enddo                 !ilm1
          endif
          enddo                 !io1
        enddo                   !ig
      enddo                     !ib1
C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_FSMBPW,procid,"fsmbpw")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endifC
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
C#endifC
C      allocate(buffer(1:3*nbas*numq), stat=ierr)
C      call MPI_ALLREDUCE(xf,buffer,3*nbas*numq,
C     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
C      if (mlog) then
C        call gettime(datim)
C        call awrit3(' fsmbpw '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' allreduce f 3*nbas*numq=%i',' ',256,lgunit(3),
C     .        procid,numprocs,3*nbas*numq)
C      endif
C      call daxpy(3*nbas*numq,1d0,buffer,1,f,1)
C      deallocate(buffer, stat=ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
C#endifC
C      deallocate(index, stat=ierr)
C      deallocate(xf, stat=ierr)
C#endif

      call tcx ('fsmbpw')


      end
