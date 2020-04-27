      subroutine drvmelp2( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx_in,nband,itp0,
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       rmelt,cmelt)
      use m_readqg, only : readqg
      use m_readeigen,only : readgeig
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
c-----------------------------------------------------
      implicit none
      real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4):: ngp1, ngp2, ngpmx_in,nqbz,nqibz,  !iq,kp,k,
     &      ngcmx ,nctot,nband,
     &      ntp0,nt0,nbloch,noccx, ngc,
     &             ngvecpB1(3,ngpmx_in),
     &             ngvecpB2(3,ngpmx_in),
     &             itp0(ntp0),
     &     nadd(3),isp  !,ngpn(nqbz)
C     integer(4):: ngveccBr(3,ngcmx)
      complex(8),allocatable::  zmelpl(:,:,:)
C     complex(8),allocatable::  geigq(:,:)
      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0)
     &         ,cmel(nbloch,noccx,ntp0)
     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)

!      complex(8) :: geigB  (ngpmx_in,nband,nqbz)
      complex(8),allocatable::geig1(:,:),geig2(:,:)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
      logical:: debug=.false.
      if(verbose()>=90) debug=.true.

      allocate( geig1(ngpmx_in,nband),geig2(ngpmx_in,nband))
c
c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) print *,' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,   ginv, qt,   ngp1, ngvecpB1) !qt is used q in practice.
      call readqg('QGpsi', q_rk,ginv, q_rkt,ngp2, ngvecpB2)
c      write(*,"('ddd1: ',i4,3f10.6,2x,3f10.6,3x,i15)")ngp1,q,qt      ,sum(abs(ngvecpB1(:,1:ngp1)))
c      write(*,"('ddd1: ',i4,3f10.6,2x,3f10.6,3x,i15)")ngp2,q_rk,q_rkt,sum(abs(ngvecpB2(:,1:ngp2)))
      if(debug) print *,'drvmelp2: goto readgeig'
      call readgeig(q,    ngpmx_in, isp, qu1, geig1)
      call readgeig(q_rk, ngpmx_in, isp, qu2, geig2)
      if(debug) print *,'drvmelp2: endof readgeig'
c      qu1=qt
c      qu2=q_rkt
c      geig1=0d0
c      geig2=0d0
c      write(*,"('ddd2: ',i4,3f10.6,2x,3f10.6,3x,d23.16,i15)")ngp1,q,qt      ,sum(abs(
c     &       geig1(1:ngp1,itp0(1:ntp0)))), ntp0  ! q1=q    ---> iq ngp1 1:ntp0 q-point
c      write(*,"('ddd2: ',i4,3f10.6,2x,3f10.6,3x,d23.16,i15)")ngp2,q_rk,q_rkt,sum(abs(
c     &       geig2(1:ngp2,1:nt0))), nt0     ! q2=q-rk ---> kp ngp2 1:nt0  occupied
cc      allocate( zmelpl(ngc,nt0,ntp0) )
c      print *,'test 888888888888888888888',ngc,nt0,ntp0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Cstop2rx 2013.08.09 kino      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(qt-qu1))>1d-10) call rx( 'drvmelp2;qu1/=qu1x')
Cstop2rx 2013.08.09 kino      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
      if(sum(abs(q_rkt-qu2))>1d-10) call rx( 'drvmelp2;qu2/=qu2x')
      if(debug) print *,' drvmelp2:end of read geig '

c     ngc  = ngcni(k)  ! k-points in IBZ
c      allocate( geigq(ngp1,ntp0) )
c      do itx = 1,ntp0
c        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c      enddo

c     print *,' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp)
      qdiff = matmul(symope,qik)  - qt + q_rkt
      ! rk    -q  +(q-rk) is not zero.
!!    <rk q-rk |q>
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) then
        print *,' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)')
     &  add,nadd,abs(add-nadd)
Cstop2rx 2013.08.09 kino        stop "drvmlp: abs(add-nadd))>1d-10"
        call rx( "drvmlp: abs(add-nadd))>1d-10")
      endif

      if(debug) print *,' goto alloc zmelpl'
      allocate( zmelpl(ngc,nt0,ntp0) )
cccccccccccccccccccccccccccccccccccccccccccccccc
c      if(.false.) then !old version
c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)), qbas,
c     o       zmelpl)
c      else             !new version from Sep-10 2000
c       call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied
c     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas,
c     o       zmelpl)
c       endif
      if(debug) print *, 'drvmelp2:goto melpln2'
ccccccccccccccccccc
c      print *, 'ssss: sum check geig1=',sum(abs(geig1(1:ngp1,itp0(1:ntp0))))
c      print *, 'ssss: sum check geig2=',ngp2,nt0,sum(abs(geig2(1:ngp2,1:nt0)))
ccccccccccccccccccccccccccccccccccc


      call melpln2( ngp1, ngvecpB1
     &           ,  ngp2, ngvecpB2
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,! q1=(shifted q) ->iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,1:nt0), nt0, ! q2=(shifted q-rk) -> kp ngp2 1:nt0  occupied
     i       shtv, matmul(symope,qik),qik, symope, qbas,
     o       zmelpl)
      if(debug) print *,'end of melpln2'

 8888 continue
c     deallocate(geig1,geig2)
      rmelt=0d0
      cmelt=0d0
      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =
     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
cc extensiton to nbloch+ngc
      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )
      deallocate(zmelpl)
      if(debug) print *,'end of drvmelp2'
      end

      subroutine melpln2
     &       ( ngp1, ngvecp1, ngp2, ngvecp2, ngc, nadd,
     &       geigq1, ntp0, ! q1=q    ---> iq 1:ntp0 q-point
     &       geigq2, nt0,  ! q2=q-rk ---> kp 1:nt0  occupied
     i       shtv,q, qi, symope, qbas,
     o       zmelp)
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
co zmelp
cr use PPOVL through getppx in m_read_ppovl
c-------------------------------------------------
      use m_read_ppovl,only:getppx, ppx,ngc2,ngvecc2
      implicit none
      integer(4):: ngp1, ngvecp1(3,ngp1), ngp2, ngvecp2(3,ngp2),
     &     itp,it, ngc,ig1,ig2,ntp0,nt0,
     &     nadd(3)   !ngvecc(3,ngc)
      complex(8):: geigq1(ngp1,ntp0), geigq1_(ntp0,ngp1),
     &          geigq2(ngp2,nt0),ge2,  !it was zmelp(ngc,nt0,ntp0); but no effect sep2009.
     &             img=(0d0,1d0) ,phas
      real(8)::  q(3),qi(3),qg(3),qbas(3,3),shtv(3)
     & ,pi=3.1415926535897932D0, symope(3,3)
      integer(4),allocatable :: ipp(:,:,:)
      complex(8),allocatable :: gg(:,:),gg_(:,:,:)
      integer(4) :: igc2,ngcx,nn(1:3),
     &       nxmin,nymin,nzmin,nxmax,nymax,nzmax
      integer(4)::verbose
c      complex(8),allocatable ::zmelp(:,:)
c      allocate(zmelp(ngc,nt0*ntp0))
      complex(8):: zmelp(ngc,nt0*ntp0)

c------------------------------
      if(verbose()>=90) print *,' melpln2: '
c- Get ppx and ngvecc2
c   ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
      call getppx(qbas,  ngcx,qi,symope) ! rppovl.f
      if(verbose()>=90) print *,' melpln2: 000111'
      if(ngc/=ngcx) then
        print *, qi, ngcx, ngc
        stop 'melpln2: ngc/= ngcx by getppx:PPOVL'
      endif
      do ig1=1,ngp1
        do itp=1,ntp0
          geigq1_(itp,ig1)=geigq1(ig1,itp)
        enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' melpln2 old sumcheck ngc2 ngvec ppx ='
c     &  ,ngc2,sum(abs(ngvecc2(1:3,1:ngc2))),sum(abs(ppx(1:ngc,1:ngc2)))
c      stop 'test end'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(verbose()>=90) print *,' melpln2: aaa'
      nxmin =minval(ngvecc2(1,1:ngc2)) ;nxmax=maxval(ngvecc2(1,1:ngc2))
      nymin =minval(ngvecc2(2,1:ngc2)) ;nymax=maxval(ngvecc2(2,1:ngc2))
      nzmin =minval(ngvecc2(3,1:ngc2)) ;nzmax=maxval(ngvecc2(3,1:ngc2))
c      write(6,"(' iqi ngc ngc2=',3i8)") iqi,ngc,ngc2
c      print *,' scheck2 ppx=',sum(ppovlinv(1:ngcx,1:ngcx))
c      print *,' scheck2 ppovls   =',sum(ppovl   (1:ngcx,1:ngc2))
c      print *,' scheck2 ngvecc2s =',sum(ngvecc2(1:3,1:ngc2))
c      print *," nxmin:nxmax,nymin:nymax,nzmin:nzmax=",
c     & nxmin,nxmax,nymin,nymax,nzmin,nzmax
      allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
      ipp = -10000000
      if(verbose()>=90) print *,' melpln2: bbb'
      do igc2 =1,ngc2
        ipp(ngvecc2(1,igc2),ngvecc2(2,igc2),ngvecc2(3,igc2))=igc2
c        nn(1:3) = ngvecc2(1:3,igc2)
c        alat =10.26
c        print *,'  nn =',nn(1:3)
c        print *,'  abs=',2*pi/alat
c     &     *sqrt(sum( (qi(1:3)+ matmul(qbas, nn(1:3)))**2 ) )
      enddo

c      inum =0
c
      if(verbose()>=90) print *,' melpln2: sht nadd=',shtv,nadd
c      do igc = 1,ngc
c        qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, ngvecc(1:3,igc)))
c        phase(igc) = exp(-img*sum(qg*shtv))
c
c      print *,' igc phase=', phase(igc)
c      enddo
c      print *," melpln2:goto ig1 ig2 loop"


      allocate(gg(ngc2,nt0*ntp0),gg_(ntp0,nt0,ngc2)) !sf
      gg_=0d0 !sf
      if(verbose()>=90) print *,' melpln2: dddd'
      do ig1 = 1,ngp1
        do ig2 = 1,ngp2
          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
          if(nn(1)<nxmin.or.nn(1)>nxmax) stop "melpln2:nn(1) out of range"
          if(nn(2)<nymin.or.nn(2)>nymax) stop "melpln2:nn(2) out of range"
          if(nn(3)<nzmin.or.nn(3)>nzmax) stop "melpln2:nn(3) out of range"
        enddo
      enddo

      do ig1 = 1,ngp1
!$OMP parallel do private(ig2, nn,igc2,qg, it, ge2 , itp ,phas) shared(gg_)
        do ig2 = 1,ngp2
          nn(1:3) = ngvecp1(1:3,ig1)-ngvecp2(1:3,ig2) - nadd
          igc2    = ipp(nn(1),nn(2),nn(3))
          qg(1:3) = 2d0*pi * (q(1:3)+ matmul(qbas, nn(1:3)) )
          phas = exp(-img*sum(qg*shtv))
          do it =1,nt0
            ge2=dconjg(geigq2(ig2,it)) *phas
            do itp=1,ntp0
              gg_(itp,it,igc2) =gg_(itp,it,igc2)
     &        +  geigq1_(itp,ig1)*ge2
            enddo
          enddo
        enddo
!$OMP end parallel do
      enddo

      do itp=1,ntp0
        do it =1,nt0
          do igc2=1,ngc2
            gg(igc2,it+(itp-1)*nt0)=gg_(itp,it,igc2)
          enddo
        enddo
      enddo
      deallocate(gg_)
C     call matm(ppx,gg,zmelp,ngc,ngc2,ntp0*nt0) !ppx= matmul(ppovlinv(ngc,ngc),ppovl(ngc,ngc2)
      call zgemm ('N','N',ngc, ntp0*nt0, ngc2, (1d0,0d0), ppx, ngc, gg, ngc2, (0d0,0d0), zmelp, ngc )
      deallocate(ppx,gg,ngvecc2,ipp)
C     <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
      do  itp = 1, ngc*nt0*ntp0
        zmelp(itp,1) = dconjg(zmelp(itp,1))
      enddo
      end
