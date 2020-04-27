      subroutine fsmbl(nbas,s_site,s_spec,s_lat,
     .  vavg,q,ndimh,nlmto,iprmb,numq,nevec,evl,evec,ewgt,f)
C- Force from smoothed hamiltonian (constant potential) and overlap
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
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  vol alat plat qlat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  hhigbl phhigb hklbl gklbl fklbl hsmbl
Ci Inputs
Ci   nbas  :size of basis
Ci   vavg  :constant potential (MT zero) to be added to h
Ci   q     :Bloch wave vector
Ci   ndimh :dimension of hamiltonian
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   numq  :number of trial Fermi levels
Ci   nevec :number of occupied eigenvectors
Ci   evl   :eigenvalues
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Co Outputs
Co   f
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Jul 08 Decouple ndimh from nlmto, for PW basis
Cu   10 Apr 02 Redimensioned eh,rsmh to accommodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   10 Sep 01 Extended to local orbitals.
Cu   23 May 00 Adapted from nfp fsm_q.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPI
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      include "mpif.h"
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
      integer nbas,ndimh,nlmto,nevec,numq,iprmb(nlmto)
      double precision evl(ndimh),f(3,nbas,numq),ewgt(numq,nevec),vavg,
     .  q(3)
      double complex evec(ndimh,ndimh)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer nlms,k0,n0,nkap0
      parameter (nlms=25, k0=1, n0=10, nkap0=4)
      integer i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,iq,is1,is2,l1,l2,ik1,ik2,
     .  ivec,m,nglob,nlm1,nlm2
      integer lh1(nkap0),lh2(nkap0),nkap1,nkap2,nlm21,nlm22,nlm11,nlm12
      integer norb1,ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),
     .  blks1(n0*nkap0),ntab1(n0*nkap0)
      integer norb2,ltab2(n0*nkap0),ktab2(n0*nkap0),offl2(n0*nkap0),
     .  blks2(n0*nkap0),ntab2(n0*nkap0)
      double precision e1(n0,nkap0),rsm1(n0,nkap0),p1(3),wt,
     .                 e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx(n0)
      double complex  s(nlms,nlms,0:k0,nkap0,nkap0)
      double complex ds(nlms,nlms,0:k0,3,nkap0,nkap0),ccc,sum
C#ifdefC MPI
C      double precision, dimension(:,:,:), allocatable :: xf
C      double precision, dimension(:),     allocatable :: buffer
C      integer, dimension(:,:), allocatable :: index
C      integer iloop,ib
C#endif

C --- Setup ---
      if (nevec <= 0) return
      call tcn ('fsmbl')

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

C --- Loop over first and second site indices ---
C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_FSMBL,procid,"fsmbl")
C#endifC
C      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
C      call dstrbp(nbas,numprocs,-1,index(0,0))
C      do  iloop = 1, index(procid,0)
C        ib1 = index(procid,iloop)
C        if (mlog) then
C          call gettime(datim)
C          call awrit4(' fsmbl '//datim//' Process %i of %i on '
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
C       Row info telling fsmbl where to poke s0 made by hhibl
        call orbl(ib1,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
        call gtbsl1(4+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)

        do  ib2 = ib1+1, nbas
          is2 = s_site(ib2)%spec
          p2 = s_site(ib2)%pos
          call uspecb(0,2,s_spec,is2,is2,lh2,rsm2,e2,nkap2)
C         Column info telling fsmbl where to poke s0 made by hhibl
          call orbl(ib2,0,nlmto,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
          call gtbsl1(4+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)

C     ... M.E. <1> and <KE> between all envelopes connecting ib1 and ib2
          do  i1 = 1, nkap1
            do  i2 = 1, nkap2
              nlm1 = (lh1(i1)+1)**2
              nlm2 = (lh2(i2)+1)**2
              if (nlm1 > nlms .or. nlm2 > nlms)
     .          call rx('fsmbl: increase nlms')
              call hhigbl(11,p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),
     .          e2(1,i2),nlm1,nlm2,1,nlms,nlms,k0,s_lat%cg,s_lat%indxcg,
     .          s_lat%jcg,s_lat%cy,s_lat,
     .          s(1,1,0,i1,i2),ds(1,1,0,1,i1,i2))
            enddo
          enddo

C     ... Loop over pairs of orbital groups, multiply bloc of gradients
          do  io2 = 1, norb2
          if (blks2(io2) /= 0) then
C           l2,ik2 = l and kaph indices, needed to locate block in s0
            l2  = ltab2(io2)
            ik2 = ktab2(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks2(io2)-1
            do  io1 = 1, norb1
            if (blks1(io1) /= 0) then
C             l1,ik1 = l and kaph indices, needed to locate block in s0
              l1  = ltab1(io1)
              ik1 = ktab1(io1)
              nlm11 = l1**2+1
              nlm12 = nlm11 + blks1(io1)-1

C         ... Loop over occupied eigenstates and x,y,z
              do  ivec = 1, nevec
              do  m = 1, 3

C           ... Loop over orbital pairs within the groups
                sum = 0d0
C               i2 = hamiltonian offset
                i2 = offl2(io2)
                do  ilm2 = nlm21, nlm22
                  i2 = i2+1
C                 i1 = orbital index in iprmb order
                  i1 = offl1(io1)
                  do  ilm1 = nlm11, nlm12
                    i1 = i1+1
                    ccc = vavg*ds(ilm1,ilm2,0,m,ik1,ik2)
     .                  -      ds(ilm1,ilm2,1,m,ik1,ik2)
     .                  - evl(ivec)*ds(ilm1,ilm2,0,m,ik1,ik2)
                    sum = sum + dconjg(evec(i1,ivec))*ccc*evec(i2,ivec)
                  enddo
                enddo
                do  iq = 1, numq
                  wt = ewgt(iq,ivec)
C#ifdefC MPI
C                  xf(m,ib1,iq) = xf(m,ib1,iq) - 2*wt*sum
C                  xf(m,ib2,iq) = xf(m,ib2,iq) + 2*wt*sum
C#else
                  f(m,ib1,iq) = f(m,ib1,iq) - 2*wt*sum
                  f(m,ib2,iq) = f(m,ib2,iq) + 2*wt*sum
C#endif
                enddo
              enddo
              enddo

            endif
            enddo
          endif
          enddo

        enddo
      enddo

C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_FSMBL,procid,"fsmbl")
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
C        call awrit3(' fsmbl '//datim//' Process %i of %i on '
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
C#endif

C      write(6,*) '---- END OF FSMBL ---'
C      do  ib = 1, nbas
C        do  iq = 1, numq
C          write(stdo,220) ib,iq,(f(m,ib,iq),m=1,3)
C  220     format(2i5,3f12.6)
C        enddo
C      enddo

      call tcx ('fsmbl')

      end
