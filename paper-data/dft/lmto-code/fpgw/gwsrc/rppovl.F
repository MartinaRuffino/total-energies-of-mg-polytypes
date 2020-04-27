      module m_read_ppovl
!! == read PPOVLGG,PPOVLG,PPOVLI ==
!!  Subroutine getppx2 returns (with common blocks in this module):
!!  nggg   :number of G vectors Gc+G2-G1  (file PPOVLGG)
!!  ggg    :<Gc G2|G1>, equivalent to <-Gc+G1-G2> (file PPOVLGG)
!!  nvggg  :vector of Gc+G2-G1 as multiples of qlat (file PPOVLGG)
!!  ngcgp  :number of Gc (?) (file PPOVLGG)
!!  nvgcgp2:vector of symope(Gc)+G2 as multiples of qlat (?) (file PPOVLGG)
!!  ngcread:number of IPWs in product basis as in: (file PPOVLGG)
!!  ngvecc :vectors of IPWs as multiples of qlat (file PPOVLGG)
!! ppovlinv:inverse of <Gc|Gc> overlap matrix (file PPOVLG.iq)
!!  Subroutine getppx2 returns (with common blocks in this module):
!! == read PPOVL ===
!!  Subroutine getppx returns (with common blocks in this module):
!!  ngc2   :number of G vectors at a give q representing eigenfunctions (file PPOVL)
!!  ngvecc2:G vectors at a given qi as multiples of qlat (file PPOVL)
!!  ppx    :ppx(1:ngc,1:ngc2) overlap between IPWs for G product basis and eigenfunction G vectors
      implicit none
      integer(4):: ngc2
      complex(8),allocatable :: ppx(:,:) ! Read from PPOVL
      integer(4),allocatable :: ngvecc2(:,:) ! Read from PPOVL
      complex(8),allocatable:: ggg(:) ! Read from PPOVLGG
      integer,allocatable:: nvggg(:,:),nvgcgp2(:,:) ! Read from PPOVLGG
      complex(8),allocatable:: ppovlinv(:,:) ! Read from PPOVLG.iq
      integer,allocatable:: ngvecc(:,:) ! Read from PPOVLG.iq
      integer:: nggg,ngcgp,ngcread
!! private
      integer(4),private:: iqix=-1, ippovl=0, ngcmx, ngc2mx, nqq, ngggmx, ngcgpmx,nqini,nqnumt
      logical,private:: ppovlclosed=.true.,init=.true.
      integer(4),allocatable,private :: ngcx_s(:),ngc2_s(:) ,ngvecc2_0_s(:,:,:)
      real(8),allocatable,private    :: qx_s(:,:)
      complex(8),allocatable,private :: ppx_s(:,:,:)
      logical,private:: debug=.false.
c-----------------------------------------------------
      contains
c-----------------------------------------------------
      subroutine getppx2(qbas,qi)
C- Return nvggg,nvgcgp2,ngvecc,  nggg,ngcgp,ngcread, ggg,ppovlinv
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovlg = 9221,ippovli = 9222,idxk
      real(8)::qbas(3,3),qi(3),qx(3)
      complex(8),allocatable :: ppxs(:,:,:)
      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:)
      real(8),allocatable :: qxs(:,:)
      logical :: KeepPPOVL,debug=.false., ippovlggooo=.true.
      integer:: verbose,io,ippovlgg=9223,ngcread2
      character*3 charnum3
!! because of get_ngpgp iqix can be not consistent among PPOVLG and PPOVLI --> use iqix=0
      if(verbose()>=100) debug=.true.
      if(ippovlggooo) then
         open(ippovlgg,file= "PPOVLGG",form='unformatted')
         rewind ippovlgg
         read(ippovlgg) nggg, ngcgp, nqq, nqini,nqnumt
         write(6,"('Readin getppx2: nggg ngcgp nqq=',3i5)") nggg, ngcgp, nqq
         allocate(nvggg(1:3,1:nggg),ggg(1:nggg),nvgcgp2(1:3,ngcgp))
         read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
         read(ippovlgg) nvggg(1:3,1:nggg)
         read(ippovlgg) ggg(1:nggg)
         close(ippovlgg)
         ippovlggooo=.false.
      endif
      do iqi = nqini,nqnumt
         open(ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
         rewind ippovlg
         read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
         if(debug) print *,'getppx2: iqi qi=',iqi,qi,qx
         ngc = ngcread
         if(sum(abs(qx-qi))>1d-10) then
           close(ippovlg)
         else
           open(ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
           rewind ippovli
           read(ippovli) qx, ngcread2 !, ngcx_s(iqi),ngc2_s(iqi)
           if(ngcread/=ngcread2) call rx('rppovl.F: inconsistent PPOVLI PPOVLg')
           if(ngc==0) exit
           if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
           allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
           read(ippovlg) ngvecc(1:3,1:ngc)
           read(ippovli) ppovlinv(1:ngc,1:ngc)
           close(ippovlg)
           close(ippovli)
           goto 1012
         endif
      enddo
      write(6,"('qi qx=',3f13.5,3x,3f13.5)") qi,qx
      call rx('getppx2: can not find given qi')
c$$$
c$$$      iqix=0
c$$$      do io=1,2 !The 1st path start from saved iqix. The second path(io=2) start from iqix=0.
c$$$        if(iqix<=0) then
c$$$           if(verbose()>=100) print *,'iqix <0 readin PPOVLG'
c$$$           if(iqix==0) close(ippovlg)
c$$$           if(iqix==0) close(ippovli)
c$$$           open(ippovlg,file= "PPOVLG",form='unformatted')
c$$$           open(ippovli,file= "PPOVLI",form='unformatted')
c$$$           iqix=0
c$$$        endif
c$$$        do iqi = iqix+1,nqq
c$$$           if(debug) print *,'getppx2: iqi qi=',iqi,qi
c$$$           read(ippovli) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
c$$$           read(ippovlg) qx, ngcread !, ngcx_s(iqi),ngc2_s(iqi)
c$$$           if(debug) print *,'getppx2: iqi qx=',iqi,qx
c$$$           ngc=ngcread
c$$$           if(sum(abs(qx-qi))<1d-10) then !hit
c$$$              if(ngc/=0) then
c$$$               if(allocated(ppovlinv)) deallocate(ppovlinv,ngvecc)
c$$$               allocate(ppovlinv(1:ngc,1:ngc),ngvecc(1:3,1:ngc))
c$$$               read(ippovlg) ngvecc(1:3,1:ngc)
c$$$               read(ippovli) ppovlinv(1:ngc,1:ngc)
c$$$              endif
c$$$              iqix=iqi+1
c$$$              goto 1012
c$$$           else
c$$$              read(ippovlg)
c$$$              read(ippovli)
c$$$           endif
c$$$        enddo
c$$$        if(iqix==0) exit
c$$$        iqix=0
c$$$      enddo
c$$$      call rx('getppx2: can not find given qi')
 1012 continue
c      print *,'eeeeeeee111111 iqix=',iqix
!! rotation of gvector rotgvec is moved to melpl2t nov2013
      end subroutine getppx2

      subroutine getppx(qbas,ngc,qi,symope)
C- Return ppx and ngvecc2
      implicit none
      integer(4):: ngc,ngcmx, ngc2mx, nqibznq0i, iqi,ngcx,ippovl = 9221,idxk
      real(8)::qbas(3,3),qi(3),symope(3,3),qx(3)
      complex(8),allocatable :: ppxs(:,:,:)
      integer(4),allocatable :: ngvecc2s(:,:,:),ngvecc2s0(:,:,:),
     &   ngcs(:), ngc2s(:),ngvecc2_0(:,:)
      real(8),allocatable :: qxs(:,:)
      logical :: KeepPPOVL
      integer:: verbose
      if(debug .and. allocated(ngvecc2_0_s)) then
        print *,' getppx :sumcheck ngvecc2_0= 3',ngc2_s(3),sum(abs(ngvecc2_0_s(1:3,1:ngc2_s(3),3)))
      endif
c--- KeepPPOVL=F to reduce memory usage.
      if(allocated(ppx)    ) deallocate(ppx)
      if(allocated(ngvecc2)) deallocate(ngvecc2)
      ngc2=999999
      if(init) then
        if(KeepPPOVL()      ) print *,' KeepPPOVL=T: m_read_ppovl keeps PPOVL in memory'
        if(.not.KeepPPOVL())  print *,' KeepPPOVL=F: m_read_ppovl reads PPOVL from disk'
        init=.false.
      endif

      if(.not.KeepPPOVL()) then
        do
          if(ppovlclosed) then
            open(ippovl,file= "PPOVL",form='unformatted')
            read(ippovl) ngcmx, ngc2mx, nqq
            ppovlclosed =.false.
c            print *,' ppovl open =',ngcmx, ngc2mx, nqq
            iqix=0
          endif
c          print *,' getppx: ppovlclosed iqix nqq=',ppovlclosed, iqix,nqq
          do iqi = iqix+1,nqq
c            print *, ' iqi=',iqi
            read(ippovl) qx, ngcx,ngc2
c            print *, ' qx ngcx ngc2=', qx, ngcx,ngc2
            if(ngcx/=0) then !Aug2005
              allocate(ngvecc2_0(1:3,1:ngc2), ppx(ngcx,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
              read(ippovl) ngvecc2_0(1:3,1:ngc2),ppx(1:ngcx,1:ngc2)
            endif
            if(sum(abs(qx-qi))<1d-10) then
              ngc=ngcx
c              if(ngc/=ngcx) then
c                print *,iqi, qx, qi, ngcx, ngc
c                stop 'getppx: ngc/= ngc from PPOVL'
c              endif
              iqix=iqi
              goto 1110
            endif
            if(ngcx/=0) deallocate(ngvecc2_0,ppx)         !sf..5June02
          enddo
Cstop2rx 2013.08.09 kino          if(iqix==0) stop ' getppx: Wrong PPVOL?'
          if(iqix==0) call rx( ' getppx: Wrong PPVOL?')
          close(ippovl)
          ppovlclosed =.true.
        enddo
      else
        if(ppovlclosed) then
          open(ippovl,file= "PPOVL",form='unformatted')
          read(ippovl) ngcmx, ngc2mx, nqq
          if(verbose()>=100) print *,' ppovl open =',ngcmx, ngc2mx, nqq
          allocate(qx_s(3,nqq), ngcx_s(nqq),ngc2_s(nqq)
     &              ,ngvecc2_0_s(3,ngc2mx,nqq),ppx_s(ngcmx,ngc2mx,nqq))
          do iqi = 1,nqq
            if(verbose()>=100)  print *, ' iqi=',iqi
            read(ippovl) qx_s(:,iqi), ngcx_s(iqi),ngc2_s(iqi)
            ngcx = ngcx_s(iqi)
            ngc2 = ngc2_s(iqi)
            if(verbose()>=100) print *,' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
cccccccccccccccccccccccccccc
c           print *, ' qx ngcx ngc2=', qx_s(:,iqi), ngcx,ngc2
cccccccccccccccccccccccccccccc
            if(ngcx/=0) then !Aug2005
C              print *,' reading ppx_s'
C              allocate(nread(1:3,ngc2), ppxread(ngcx,ngc2))
C              read(ippovl) nread, ppxread
C              ngvecc2_0_s(1:3,1:ngc2,iqi)=nread
C              ppx_s(1:ngcx,1:ngc2,iqi)=ppxread
C              deallocate(nread,ppxread)
              print "('ngcmx,ngcx, ngc2mx,ngc2, nqq,iqi: ',6(x,i0))", ngcmx,ngcx, ngc2mx,ngc2, nqq,iqi
              read(ippovl) ngvecc2_0_s(1:3,1:ngc2,iqi),ppx_s(1:ngcx,1:ngc2,iqi)
            endif
          enddo
          ppovlclosed =.false.
          if(verbose()>=100) print *,' end of readin PPOVL'
        endif
        if(verbose()>=100) print *,' goto iqi search'
ccccccccccccccccccccccccccccccccc
c           print *,' qi=',qi,nqq
ccccccccccccccccccccccccccccccccc
        iqi = idxk(qi,qx_s,nqq)
        ngc2= ngc2_s(iqi)
        if(verbose()>=100) print *,' iqi ngc2',iqi,ngc2
        ngc=ngcx_s(iqi)
        allocate(ngvecc2_0(3,ngc2), ppx(ngc,ngc2))  !ppx(ngc2+1 is for bug of zgemm in matm when it is combined with alloc.
        ngvecc2_0(1:3,1:ngc2)= ngvecc2_0_s(1:3,1:ngc2,iqi)
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc
        ppx(1:ngc,1:ngc2)   = ppx_s(1:ngc,1:ngc2,iqi)
c          if(ngc/=ngcx_s(iqi)) then
c            print *,iqi, qx, qi, ngcx, ngc
c            stop 'getppx: xxx ngc/= ngc from PPOVL'
c          endif
      endif
c---------------------------------------------
 1110 continue
      allocate(ngvecc2(1:3,1:ngc2),ngc2s(1))
      ngc2s(1) = ngc2
      call rotgvec(symope, 1,
c     i    ngc2mx, ngc2s,qbas, ngvecc2_0, !alpha linux caused segmentation error
     i    ngc2, ngc2s,qbas, ngvecc2_0,
     o    ngvecc2)
      if(allocated(ngc2s)    ) deallocate(ngc2s)
      if(allocated(ngvecc2_0)) deallocate(ngvecc2_0)
      end subroutine
      end module
