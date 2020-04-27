!! read BZDATA
      module m_read_bzdata
c--------readin data set when you call read_BZDATA ---------------
!! The following are output when you call read_BZDATA()
      implicit none
      integer :: n1,n2,n3,ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,itet
      integer,allocatable ::
     .  idtetf(:,:),ib1bz(:),idteti(:,:),
     .  nstar(:),irk(:,:),nstbz(:) !,index_qbz(:,:,:)
      real(8):: qbas(3,3),ginv(3,3),dq_(3)
      real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:)
      logical:: done_read_bzdata=.false.

c-----------------------------------------------------------------
      contains

      subroutine read_BZDATA()
c No input except BZDATA file.
c- After you call this, you can access Brillowin Zone datas above ----
!ML  MPI: This routine reads on rank 0 and broadcasts to the other processes
      use mpi_mod, rank => mpi_rank

      implicit none
      integer :: intq(3),iqbz,ifbz,n,verbose,i,narray(10)
      real(8) :: qout(3)
      ifbz = 6661

      if(verbose()>40) write(6,'("[",I4,"] ### readin BZDATA ###")') rank

      if (rank == 0) then
        open (ifbz, file='BZDATA')
!        read(ifbz,*)  nqbz,nqibz,nqbzw,ntetf,nteti,ngrp !,nqibz_r
!        read(ifbz,*)  n1,n2,n3
         read(ifbz,*) (narray(i),i=1,6)
         read(ifbz,*) (narray(i),i=8,10)
      endif
#ifdef USE_MPI
      call MPI_Bcast(narray,10,MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
#endif
      nqbz    = narray(1)
      nqibz   = narray(2)
      nqbzw   = narray(3)
      ntetf   = narray(4)
      nteti   = narray(5)
      ngrp    = narray(6)
!     nqibz_r = narray(7)
      n1 = narray(8)
      n2 = narray(9)
      n3 = narray(10)
      if(rank == 0) write(*,321) ngrp, nqibz, nteti, nqbz, ntetf
  321   format(/' read_bzdata:',i3,' symops'/
     .    i6,',',i6,' qp, tetrahedra in the irreducible BZ'/
     .    i6,',',i6,' qp, tetrahedra in the full BZ')

      allocate(qbz(3,nqbz),wbz(nqbz))
      allocate(qibz(3,nqbz),wibz(nqbz),nstbz(nqbz))
      allocate(nstar(nqibz), irk(nqibz,ngrp))
      if(ntetf>0) then
        allocate( idtetf(0:3,ntetf), ib1bz(nqbzw), qbzw(3,nqbzw) )
      endif
      if(nteti>0) then
        allocate( idteti(0:4,6*nqbz))
      endif
C     allocate( qibz_r(3,nqibz_r) )

      if (rank == 0) then
        call rwbzdata(ifbz,1,
     &    ngrp,qbas,ginv,
     i    qbz,wbz,nstbz,nqbz,
     i    qibz,wibz,nstar,irk,nqibz,
     i    idtetf,ntetf,qbzw,ib1bz,nqbzw,
     i    idteti,nteti,dq_)
      close(ifbz)
      endif
#ifdef USE_MPI
      call MPI_Bcast(qbas,  9,MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ginv,  9,MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
!     call MPI_Bcast(qbasmc,9,MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)

      call MPI_Bcast(qibz,  3*nqbz, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(wibz,  nqbz,   MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nstar, nqibz,   MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
      call MPI_Bcast(irk,   nqibz*ngrp, MPI_INTEGER,0, MPI_COMM_WORLD,ierror)

!     call MPI_Bcast(nqibz_r, 1, MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
!     call MPI_Bcast(qibz_r, 3*nqibz_r, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)

      call MPI_Bcast(qbz, 3*nqbz, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(wbz,   nqbz, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nstbz, nqbz, MPI_INTEGER,0, MPI_COMM_WORLD,ierror)

      if(ntetf>0) then
         call MPI_Bcast(idtetf,4*ntetf, MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
         call MPI_Bcast(ib1bz,  nqbzw,  MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
         call MPI_Bcast(qbzw, 3*nqbzw,  MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      endif

      if(nteti>0) then
         call MPI_Bcast(idteti, 5*6*nqbz, MPI_INTEGER,0, MPI_COMM_WORLD,ierror)
      endif

      call MPI_Bcast(dq_, 3, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
#endif

C ... Add index to specify qbz. index_qbz is used to find index in qbz.
C This algolism here can make index_qbz for any qbz(1:3,nqbz).
c With fbz2, you can get index iq for given q(1:3) so that
c            q(1:3)= qbz(1:3,iq) + some G vector.
c
c To find iqbz for given q,
c      call rangedq(matmul(ginv,q), qout)
c      intq =  qout*n +1
c      iqbz= index_qbz(intq(1),intq(2),intq(3))
c See fbz2.

c         n = 0
c         allocate(index_qbz(1,1,1)) !dummy
c 1120    continue
c         n = n+1
c         print *,' =========================== n=',n
c         deallocate(index_qbz)
c         allocate(index_qbz(n,n,n))
c         index_qbz = -9999
c         do iqbz = 1,nqbz
c           call rangedq(matmul(ginv,qbz(:,iqbz)), qout)
c           intq =  qout*n +1
cc           print *,' qbz=', qbz(:,iqbz) !,matmul(ginv,qbz(:,iqbz))
cc           print *,' qout='
cc           print *, qout
cc           print *, ' intq=', intq
c           if(verbose()>=100) write(6,"(' qbz=',3f10.5,'  index=',3i3,'   qout=',3f10.5)")
c     &       qbz(:,iqbz),intq,qout
c           call checkrange(intq(1),1,n) !sanity checks
c           call checkrange(intq(2),1,n)
c           call checkrange(intq(3),1,n)
c           if(index_qbz(intq(1),intq(2),intq(3))/=-9999) then
c              if(verbose()>=100) print *,'failed indexing with this n. Try to enlarged n'
c              goto 1120
c           endif
c           index_qbz(intq(1),intq(2),intq(3)) = iqbz
c         enddo
c         n_index_qbz=n
c         write(6,*) " O.K. index_qbz is generated; n=", n

      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
        print *, 'sum (wibz)=', sum(wibz(1:nqibz))
        call rx( 'read_BZDATA  sum (wibz) is not 2.')
      endif
      done_read_bzdata=.true.

      end subroutine read_BZDATA
      end module

      subroutine rwbzdata(ifbz,job,
     .  ngrp,qbas,ginv,
     .  qbz,wbz,nstbz,nqbz,
     .  qibz,wibz,nstar,irk,nqibz,
     .  idtetf,ntetf,qbzw,ib1bz,nqbzw,
     .  idteti,nteti,dq_)
C- Read/Write BZ mesh data required for GW
C----------------------------------------------------------------
      implicit none
      integer:: nqbz,ntetf,nteti,nqbzw,iqbz,ifbz
     & ,nqibz,iqibz,itet,ngrp,job
      real(8):: qbz(3,nqbz),wbz(nqbz),qibz(3,nqibz),wibz(nqibz)
     &         ,qbzw(3,nqbzw),qbas(3,3),ginv(3,3),dq_(3)
      integer:: idtetf(0:3,ntetf),ib1bz(nqbzw),idteti(0:4,nteti)
     &       ,irk(nqibz,ngrp),nstar(nqibz),nstbz(nqbz)
      character(len=16) :: irk_fmt

      write(irk_fmt, "('(',i0,'(x,i8))')") ngrp

      if(job<0) write(ifbz,"(3d24.16)") qbas,ginv
      if(job>0) read (ifbz,"(3d24.16)") qbas,ginv

      do iqibz = 1,nqibz
        if(job<0) then
          write(ifbz,"(4d24.16,i9)")
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          write(ifbz,irk_fmt) irk(iqibz,1:ngrp)
        else
          read(ifbz,"(4d24.16,i9)")
     &    qibz(1:3,iqibz),wibz(iqibz),nstar(iqibz)
          read(ifbz,irk_fmt) irk(iqibz,1:ngrp)
        endif
      enddo

      do iqbz = 1,nqbz
        if(job<0) then
          write(ifbz,"(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
        else
          read(ifbz, "(4d24.16,i10)") qbz(1:3,iqbz),wbz(iqbz),nstbz(iqbz)
        endif
      enddo

      if (ntetf>0) then
        if(job<0) then
          write(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          write(ifbz,"(i9,3d24.16)") (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        else
          read(ifbz,"(4i10)") (idtetf(0:3,itet),itet=1,ntetf)
          read(ifbz,"(i9,3d24.16)") (ib1bz(iqbz), qbzw(1:3,iqbz),iqbz=1,nqbzw)
        endif
      endif

      if(nteti>0) then
        if(job<0) then
          write(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        else
          read(ifbz,"(5i10)") (idteti(0:4,itet),itet=1,nteti)
        endif
      endif

      if(job<0) write(ifbz,"(3d24.16,' !dq_')") dq_
      if(job>0) read (ifbz,"(3d24.16)") dq_

      end subroutine rwbzdata

      subroutine checkrange(intq,n1,n2)
      if(intq<n1 .or. intq>n2) then
        print *,'checkrange: intq n1 n2= ',intq,n1,n2
        call rx( 'checkrange: stop ')
      endif
      end
