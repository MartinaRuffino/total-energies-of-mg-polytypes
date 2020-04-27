      module m_readq0p
      real(8),allocatable,protected:: wqt(:), wgt0(:,:),q0i(:,:) !,nx(:,:),nblocha(:)
      integer,protected:: nq0i,nq0iadd
      integer,protected,allocatable:: ixyz(:)

      contains
      subroutine readq0p()
      use mpi_mod, rank => mpi_rank
      implicit none
      integer:: neps,ifiq0p,ifile_handle,i,nq0ix,iq0pin
      logical:: debug=.false.
c      write(6,*) 'reading QOP'
      ifiq0p=ifile_handle()

      if (rank == 0) then
        open (ifiq0p,file='Q0P')
        nq0iadd = 0 ! for now
        read (ifiq0p,*) nq0i,iq0pin,nq0iadd
        allocate( wqt(1:nq0i),q0i(1:3,1:nq0i+nq0iadd),ixyz(nq0i+nq0iadd) )
        do i=1,nq0i+nq0iadd
          read (ifiq0p, * ) wqt(i),q0i(1:3,i),ixyz(i)
        enddo
      endif
#ifdef USE_MPI
      call MPI_Bcast(nq0i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(nq0iadd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      if (rank /= 0) allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      call MPI_Bcast(wqt,nq0i,MPI_Double_Precision,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(q0i,3*nq0i,MPI_Double_Precision,0,MPI_COMM_WORLD,ierror)
#endif

      nq0ix = nq0i
      do  i=1,nq0i
        if (wqt(i) == 0d0) then
          nq0ix = i-1
          exit
        endif
      enddo
      neps = nq0i - nq0ix       ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      if (rank == 0) then
        write(6,*) ' num of zero weight q0p=',neps
        write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i+nq0iadd)
        close(ifiq0p)
      endif
      end subroutine
      end module

      subroutine readppovl0(q,ngc,ppovl)
      implicit none

      integer, intent(in) :: ngc
      complex(8), intent(out) :: ppovl(ngc,ngc)
      real(8), intent(in) :: q(3)
      integer:: ngc_r,ippovl0,ifile_handle
      real(8):: qx(3),tolq=1d-6
C      integer:: ngc_r,ippovl0,ixx,ngc
C      complex(8):: ppovl(ngc,ngc)
C      real(8):: q(3),qx(3)
      ippovl0=2301
      open(ippovl0,file='PPOVL0',form='unformatted')
      do
        read(ippovl0) qx,ngc_r
        if(sum(abs(qx-q))<tolq) then
          if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
          read(ippovl0) ppovl
          exit
        endif
      enddo
      close(ippovl0)
      end

      module m_readQG
C- Return QGcou and QGpsi
c
c     Notes about parallelization:
c
c     The routines seem  OK for MPI at the moment (as long as
c     the files and arrays are not getting too big).
c
c     Possible improvements:
c
c     1) Read only on master node and distribute arrays
c
c     2) Distribute only the data needed by a specific iq-group
c
c     All these are private.

      implicit none
      real(8),allocatable,private,target:: qc(:,:),qp(:,:)
      logical,private:: init(2)=.true.
      real(8),private:: QpGcut_cou, QpGcut_psi
      integer,private,target::   nqnumc,nqnump,ngcmx,ngpmx
      integer,allocatable,private:: ngvecp(:,:,:),ngp(:),ngvecc(:,:,:),ngc(:)
!! --------------------------------------------------------
!! See readeigen.F
c      module m_iqindx_qtt
c      use m_hamindex, only: qtt,nqtt
      integer,pointer,private::nqtt
      real(8),pointer,private::qtt(:,:)
      real(8),private:: epsd=1d-7
      integer,private,pointer:: nkey(:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      integer,target,private :: nkeyp(3),nkeyc(3)
      integer,target,allocatable,private:: keyp(:,:),kk1p(:),kk2p(:),kk3p(:),iqkkkp(:,:,:)
      integer,target,allocatable,private:: keyc(:,:),kk1c(:),kk2c(:),kk3c(:),iqkkkc(:,:,:)
      real(8),private:: ginv_(3,3)
      contains
c----------------------------------
      subroutine readngmx(key,ngmx)
c- get ngcmx or mgpmx
c
c -- MPI remark: needs to be called on all nodes !
      use mpi_mod, rank => mpi_rank

      implicit none

      integer:: ngmx,ifiqg=4052
      character*(*) key
      logical :: wrong_key=.false.

      if(rank == 0) then
        if (key=='QGpsi') then
          open(ifiqg, file='QGpsi',form='unformatted')
          read(ifiqg) nqnump, ngpmx, QpGcut_psi
          ngmx = ngpmx
        elseif(key=='QGcou') then
          open(ifiqg, file='QGcou',form='unformatted')
          read(ifiqg) nqnumc, ngcmx, QpGcut_cou
          ngmx = ngcmx
        else
          wrong_key = .true.
        endif
        close(ifiqg)
      endif
#ifdef USE_MPI
      if (key=='QGpsi') then
         call MPI_Bcast(nqnump,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierror)
         call MPI_Bcast(ngpmx,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierror)
      elseif (key=='QGcou') then
         call MPI_Bcast(nqnumc,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierror)
         call MPI_Bcast(ngcmx,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierror)
      endif
      call MPI_Bcast(ngmx,1,MPI_INTEGER,0, MPI_COMM_WORLD, ierror)
#endif
      if(wrong_key) call rx( "readngmx: key is not QGpsi QGcou")

      end subroutine

      subroutine readqg(key,qin,ginv,  qu,ngv,ngvec)
C- Get ngv and ngvec(3,ngv) for given qin(3)
Cr key should be one of 'QGcou' or 'QGpsi'
      implicit none
C     I/O arguments
      character*(*), intent(in) :: key
      real(8), intent(in) :: qin(3),ginv(3,3)
      real(8), intent(out) :: qu(3)
      integer, intent(out) :: ngv, ngvec(3,*)
C     local variables
      integer ::  ifi, iq,verbose

      if (key=='QGpsi') then
        ifi=1
        if(verbose()>=80) write (6,"(' readqg psi: qin=',3f8.3,i5)") qin
      elseif (key=='QGcou') then
        ifi=2
        if(verbose()>=80) write (6,"(' readqg cou: qin=',3f8.3,i5)") qin
      else
        call rxs( "readqg: key not recognized",key)
      endif
      if (init(ifi)) then
        call init_readqg(ifi,ginv)
        init(ifi)=.false.
      endif
      call iqindx2qg(qin,ifi, iq,qu)
      if (ifi==1) then
        ngv = ngp(iq)
        ngvec(1:3,1:ngv) = ngvecp(1:3,1:ngv,iq)
        return
      elseif (ifi==2) then
        ngv = ngc(iq)
        ngvec(1:3,1:ngv) = ngvecc(1:3,1:ngv,iq)
        return
      endif
      call rx("readqg: can not find QGpsi or QPcou for given q")
      end subroutine readqg
      subroutine readqg0(key,qin,ginv,  qu,ngv)
C- Get the number of G vectors for either eigenfunction basis or mixed basis
Cr key should be one of 'QGcou' or 'QGpsi'
      implicit none
      character*(*), intent(in) :: key
      integer:: ngv, ifi, iq,verbose
      real(8):: qin(3),qu(3),ginv(3,3)
      if    (key=='QGpsi') then
        ifi=1
        if(verbose()>=80) write (6,"('readqg0 psi: qin=',3f8.3,i5)") qin
      elseif(key=='QGcou') then
        ifi=2
        if(verbose()>=80) write (6,"('readqg0 cou: qin=',3f8.3,i5)") qin
      else
        call rxs( "readqg: key not recognized",key)
      endif
      if(init(ifi)) then
        call init_readqg(ifi,ginv)
        init(ifi)=.false.
      endif

      call iqindx2qg(qin,ifi, iq,qu)
      if(ifi==1) then
c        call iqindx2(qin,ginv, qp,nqnump, iq,qu)
        ngv  = ngp(iq)
        if(verbose()>=80) print *,'ngp=',ngv
      elseif(ifi==2) then
c        call iqindx2(qin,ginv, qc,nqnumc, iq,qu)
        ngv  = ngc(iq)
        if(verbose()>=80) print *,'ngc=',ngv
      endif
      return
      call rx( "readqg0: can not find QGpsi or QPcou for given q")
      end subroutine

c--------------------------
      subroutine init_readqg(ifi,ginv)
C- Initialization for reading from files QGpsi or QGcou.
      implicit none
      integer:: ifi,ifiqg,iq,verbose
      real(8)::qq(3)
      real(8):: ginv(3,3)
      real(8),allocatable:: qxx(:,:)
      integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,ik
      integer,allocatable:: ieord(:),key(:,:)
      ginv_=ginv
      print *,' init_readqg ifi=',ifi
      ifiqg=4052
      if(ifi==1) then
        open(ifiqg, file='QGpsi',form='unformatted')
        read(ifiqg) nqnump, ngpmx, QpGcut_psi
        if(verbose()>49)
     &     write (6,"('init_readqg ngnumc ngcmx QpGcut_psi=',2i5,f8.3)")
     &     nqnump, ngpmx, QpGcut_psi
        allocate(ngvecp(3,ngpmx,nqnump),qp(3,nqnump),ngp(nqnump))
        do iq=1, nqnump
          read (ifiqg) qp(1:3,iq), ngp(iq)
          read (ifiqg) ngvecp(1:3,1:ngp(iq),iq)
          if(verbose()>40)
     &     write (6,"('init_readqg psi qp ngp =',3f8.3,i5)") qp(1:3,iq),ngp(iq)
        enddo
      elseif(ifi==2) then
        open(ifiqg, file='QGcou',form='unformatted')
        read(ifiqg) nqnumc, ngcmx, QpGcut_cou
c         write (6,"('init_readqg ngnumc ngcmx QpGcut_cou=',2i5,f8.3)")
c     &     nqnumc, ngcmx, QpGcut_cou
        allocate(ngvecc(3,ngcmx,nqnumc),qc(3,nqnumc),ngc(nqnumc))
        do iq=1, nqnumc
          read(ifiqg) qc(1:3,iq), ngc(iq)
          if(verbose()>40)  write (6,"('init_readqg cou  qc ngc =',3f8.3,i5)") qc(1:3,iq), ngc(iq)
          read (ifiqg) ngvecc(1:3,1:ngc(iq),iq)
        enddo
      endif
      close(ifiqg)

!! === mapping of qtt ===
!! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
!!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
!!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
      if(ifi==1) then
         nqtt => nqnump
         qtt  => qp
         nkey => nkeyp
      elseif(ifi==2) then
         nqtt => nqnumc
         qtt  => qc
         nkey => nkeyc
      endif
!! followings are the same as codes in readeigen.F
      allocate(ieord(nqtt))
      allocate(key(3,0:nqtt),qxx(3,nqtt))
      key(:,0)=0 !dummy
      key=-99999
      do iq=1,nqtt
         call rangedq(matmul(ginv_,qtt(:,iq)), qxx(:,iq))
      enddo
!! get key and nkey for each ix.
      do ix =1,3
         call sortea(qxx(ix,:),ieord,nqtt,isig)
         ik=0
         do i=1,nqtt
            kkk=(qxx(ix,ieord(i))+0.5*epsd)/epsd  !kkk is digitized by 1/epsd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
c               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
               call rx( 'iqindx: bug not sorted well')
            endif
         enddo
         nkey(ix)=ik
      enddo
      deallocate(ieord)
!!  key is reallocated. inverse mattping, iqkkk
      if(ifi==1) then
        allocate( kk1p(nkey(1)),kk2p(nkey(2)),kk3p(nkey(3)) )
        allocate( iqkkkp(nkey(1),nkey(2),nkey(3)) )
         iqkkk => iqkkkp
         kk1 =>kk1p
         kk2 =>kk2p
         kk3 =>kk3p
      elseif(ifi==2) then
        allocate( kk1c(nkey(1)),kk2c(nkey(2)),kk3c(nkey(3)) )
        allocate( iqkkkc(nkey(1),nkey(2),nkey(3)) )
         iqkkk => iqkkkc
         kk1 =>kk1c
         kk2 =>kk2c
         kk3 =>kk3c
      endif

      kk1(:) = key(1,1:nkey(1))
      kk2(:) = key(2,1:nkey(2))
      kk3(:) = key(3,1:nkey(3))
      deallocate(key)
c      print *,' ifi init_qqq nqtt=',ifi,nqtt
c      print *,'kkk3=',kkk3
c      print *,'nkey=',nkey
c      print *,'kk1=',kk1
      do i=1,nqtt
         kkk3= (qxx(:,i)+0.5*epsd)/epsd !kkk is digitized by 1/epsd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
c      print *,' ik1,ik2,ik3 i=',ik1,ik2,ik3,i
      enddo
      deallocate(qxx)
      end subroutine init_readqg
!! ---
      subroutine tabkk(kkin, kktable,n, nout)
      integer:: nout,n, kkin, kktable(n),i,mm,i1,i2
      i1=1
      i2=n
      if(kkin==kktable(1)) then
         nout=1
         return
      elseif(kkin==kktable(n)) then
         nout=n
         return
      endif
      do i=1,n
         mm=(i1+i2)/2
         if(kkin==kktable(mm)) then
            nout=mm
            return
         elseif(kkin>kktable(mm)) then
            i1=mm
         else
            i2=mm
         endif
      enddo
      print *, i1,i2,kkin
      print *, kktable(i1),kktable(i2)
      call rx( 'takk: error')
      end subroutine

c$$$c--- release to save memory area.
c$$$      subroutine releaseqg_notusednow(key)
c$$$      implicit none
c$$$      character*(*) key
c$$$      integer:: ifi
c$$$      if    (key=='QGpsi') then
c$$$        ifi=1
c$$$        deallocate(qp,ngvecp)
c$$$      elseif(key=='QGcou') then
c$$$        ifi=2
c$$$        deallocate(qc,ngvecc)
c$$$      else
c$$$        stop "releaseqg: in readQGcou"
c$$$      endif
c$$$      init(ifi)=.false.
c$$$      end subroutine
!!---------------------------------------------------------
      subroutine iqindx2qg(q,ifi, iqindx,qu)
!!=== Find index as q=qq(:,iq) with modulo of primitive vector.===
!! ginv is the inverse of plat (premitive translation vector).
!! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
      implicit none
      integer:: i_out, iq,iqx,iqindx ,kkk3(3),ik1,ik2,ik3,ifi
      real(8):: q(3),qx(3),qu(3),qzz(3)
      logical::debug=.false.
      if(ifi==1) then
c         nqtt => nqnump
         qtt  => qp
         nkey => nkeyp
         iqkkk => iqkkkp
         kk1 =>kk1p
         kk2 =>kk2p
         kk3 =>kk3p
      elseif(ifi==2) then
c         nqtt => nqnumc
         qtt  => qc
         nkey => nkeyc
         iqkkk => iqkkkc
         kk1 =>kk1c
         kk2 =>kk2c
         kk3 =>kk3c
      endif
      if(debug) write(*,"(' iqindx2_: q=',3f20.15)") q
      call rangedq(matmul(ginv_,q), qzz)
      kkk3 = (qzz+0.5*epsd)/epsd
c      print *,'kkk3=',kkk3
c      print *,'kk1,nkey1',kk1,nkey(1)
c      print *,'kk2,nkey2',kk2,nkey(2)
c      print *,'kk3,nkey3',kk3,nkey(3)
      call tabkk(kkk3(1), kk1,nkey(1), ik1)
      call tabkk(kkk3(2), kk2,nkey(2), ik2)
      call tabkk(kkk3(3), kk3,nkey(3), ik3)
c      print *,' ik1ik2ik3=',ik1,ik2,ik3
      iqindx = iqkkk(ik1,ik2,ik3)
c      print *,'iqindx=',iqindx
      qu = qtt(:,iqindx)
      end subroutine

!! === this is  taken from sorta.F ===
      subroutine sortea(ea,ieaord,n,isig)
c mini-sort routine.
      integer:: n,isig,ix,itmp,iwirk,i,j
      real(8)::        ea(n)
      integer:: ieaord(n)
c sorting of ea
      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
#ifdef EXPAND_ISWAP
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp

#else
            call iswap (ieaord(i-1),ieaord(i))
#endif
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      end subroutine
#ifndef EXPAND_ISWAP
      subroutine iswap(i,j)
      integer iwork, i, j
      iwork = j
      j = i
      i = iwork
      end subroutine
#endif
c---------------------------------------------------------
      end module m_readQG

      subroutine getgvec(key,qin,ginv,  qu,ngv,ngvec)
C- Get the G vectors for either eigenfunction basis or mixed basis
Cr key should be one of 'QGcou' or 'QGpsi'
Cr This routine performs the same function as readqg,  but no interface to m_readqg is required
      use m_readqg,only : readqg
      implicit none
C     I/O arguments
      integer, INTENT(OUT) :: ngv,ngvec(3,*)
      real(8), INTENT(IN)  :: qin(3), ginv(3,3)
      real(8), INTENT(OUT) :: qu(3)
      character*(*), INTENT(IN) ::  key

      call readqg(key,qin,ginv,  qu,ngv,ngvec)
      end

      subroutine getngvec(key,qin,ginv,  qu,ngv)
C- Get the number of G vectors for either eigenfunction basis or mixed basis
Cr key should be one of 'QGcou' or 'QGpsi'
Cr This routine performs the same function as readqg0,  but no interface to m_readqg is required
      use m_readqg,only : readqg0
      implicit none
C     I/O arguments
      integer, INTENT(OUT) :: ngv
      real(8), INTENT(IN)  :: qin(3), ginv(3,3)
      real(8), INTENT(OUT) :: qu(3)
      character*(*), INTENT(IN) ::  key

      call readqg0(key,qin,ginv,  qu,ngv)
      end


