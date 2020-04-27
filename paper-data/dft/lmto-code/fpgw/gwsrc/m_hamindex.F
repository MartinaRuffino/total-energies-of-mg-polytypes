      module m_hamindex
!! == this is in lm7K/subs/m_hamindex.F and in fpgw/gwsrc/m_hamindex.F ==
!! We will need to unify make system and source code in fpgw and lmf.
C  See comments at the start of readeigen.F
      implicit none
      integer,parameter:: null=-999999
      integer:: ngrp=null, lxx=null, kxx=null,norbmto=null,norbtx=null
!! norbtx is given in gwsrc/readeigen.F init_readeigen2
      integer:: imx=null,nbas,ndimham=null,nat
      integer:: nqtt,nqi, nqnum,ngpmx
      integer,allocatable:: ltab(:),ktab(:),offl(:),ispec(:), iclasst(:),offlrev(:,:,:),ibastab(:),ipbi(:)
      integer,allocatable:: iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:) !,ngvecp(:,:,:),ngvecprev(:,:,:,:)
      real(8),allocatable:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:)
      real(8):: plat(3,3),qlat(3,3)
      real(8),allocatable:: qtt(:,:),qtti(:,:)
      integer,allocatable:: igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
      integer::  napwmx=null,lxxa=null
      logical,private:: debug=.false.
      contains

      integer function getikt(qin)
      integer::i
      real(8):: qin(3)
c      if(debug) print *,'nkt=',nkt
      do i=1, nqnum !*2 !nkt
        if(debug) print *,i,qin, qq(:,i)
        if(sum (abs(qin-qq(:,i)))<1d-8) then
          getikt=i
          return
        endif
      enddo
      print *,' getikt: xxx error nqnum qin=',nqnum,qin
      do i=1, nqnum !*2 !nkt
        write(*,"('i qq=',i3,3f11.5)")i, qq(:,i)
      enddo
      call rx( ' getikt can not find ikt for given q')
      end function

      subroutine writehamindex()
      integer(4):: ifi
      logical::pmton
!! == write info for wave rotation. feb2012takao ==
      ifi=1789
      open(ifi,file='HAMindex',form='unformatted')
      write(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto
      write(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
      write(ifi)lxxa
      write(ifi)dlmm
      write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
      write(ifi)qq !,ngvecp,ngvecprev
      write(ifi)plat,qlat,napwmx
      if(napwmx/=0) then        !for APW rotation used in rotwvigg
        write(ifi) igv2,napwk,igv2rev
      endif
      write(ifi) nat, ipbi(1:nat)  ! Appended Aug 17
      close(ifi)
      end subroutine writehamindex

      subroutine readhamindex()
C- read info for wave rotation. feb2012takao
      use mpi_mod, rank => mpi_rank

      integer(4):: ifi,nkt
      logical::pmton
      ifi=1789

C     Aug 2017 MvS added nat, ipbi
      if(rank == 0) then
        open(ifi,file='HAMindex',form='unformatted')
        read(ifi)ngrp,nbas,kxx,lxx,nqtt,nqi,nqnum,imx,ngpmx,norbmto,nat
      endif
#ifdef USE_MPI
      call MPI_Bcast(ngrp,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nbas,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(kxx,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(lxx,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nqtt,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nqi,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nqnum,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(imx,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ngpmx,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(norbmto, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nat,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      allocate(symops(3,3,ngrp),ag(3,ngrp),qtt(3,nqtt),qtti(3,nqi))
      allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
      allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))

      if(rank == 0) then
        read(ifi)symops,ag,invgx,miat,tiat,shtvg,qtt,qtti,iqmap,igmap,iqimap
        read(ifi) lxxa
      endif
#ifdef USE_MPI
      call MPI_Bcast(symops, 9*ngrp, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ag,     3*ngrp, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(invgx, ngrp, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(miat, nbas*ngrp, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(tiat,3*nbas*ngrp, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(shtvg,  3*ngrp, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(qtt,    3*nqtt, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(qtti,   3*nqi,  MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(iqmap,  nqtt, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(igmap,  nqtt, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(iqimap, nqtt, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(lxxa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
      allocate( offlrev(nbas,0:lxx,kxx))
      allocate( dlmm(-lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))

      if(rank == 0) then
        read(ifi) dlmm
        read(ifi)ibastab,ltab,ktab,offl,offlrev
      endif

#ifdef USE_MPI
      call MPI_BCast(dlmm, (2*lxxa+1)**2 * (lxxa+1) * ngrp,
     &               MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(ibastab, norbmto, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(ltab,    norbmto, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(ktab,    norbmto, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(offl,    norbmto, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(offlrev, nbas*(lxx+1)*kxx,  MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

c      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
c      allocate( ngvecp(3,ngpmx,nqnum) )
      allocate( qq(3,nqnum)) !this was qq(3,nqnum*2) until Aug2012 when shorbz had been used.

      if(rank == 0) then
        read(ifi)qq             !,ngvecp,ngvecprev
        read(ifi)plat,qlat,napwmx
      endif
#ifdef USE_MPI
      call MPI_BCast(qq, 3*nqnum, MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(plat,   9,   MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(qlat,   9,   MPI_Double_Precision, 0, MPI_COMM_WORLD, ierror)
      call MPI_BCast(napwmx, 1,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      if(napwmx/=0)then !for APW rotation used in rotwvigg
        nkt=nqtt
        allocate( igv2(3,napwmx,nkt) )
        allocate( napwk(nkt))
        allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
        if(rank == 0) read(ifi) igv2,napwk,igv2rev
#ifdef USE_MPI
        call MPI_BCast(igv2, 3*napwmx*nkt, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        call MPI_BCast(napwk, nkt,         MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        call MPI_BCast(igv2rev, (2*imx+1)**3 *nkt, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif
      endif

      allocate(ipbi(nat))
      if(rank == 0) read(ifi) ipbi
#ifdef USE_MPI
      call MPI_BCast(ipbi, nat,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
#endif

      if(rank == 0) close(ifi)

      end subroutine readhamindex
      end module


