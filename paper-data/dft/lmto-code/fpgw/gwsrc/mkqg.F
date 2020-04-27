      subroutine mkQG2(alat,plat,symops,ngrp,nnn,iq0pin,
     &     QpGcut_psi, QpGcut_Cou, ifiqg, ifiqgc,gammacellctrl,lnq0iadd)
C- Make required q and G in the expansion of GW
!!     |q+G| < QpGcut_psi for eigenfunction psi.
!!     |q+G| < QpGcut_Cou for coulomb interaction
!!
!! OUTPUT
!!     file handle= ifiqg,  which contains q and G points for eigenfunction psi. --> QGpsi
!!     file handle= ifiqgc, which contains q and G points for Coulomb            --> QGcou
!!
!!     QGpsi(ifiqg), QGcou(ifiqgc), Q0P are written.
!!     See the end of console output.
! Not used: qbzwm
Cu Updates
Cu   27 Apr 19 Adapted Takao's new mkQG2
Cu   04 May 16 Minor cleanup
!! ---------------------------------------------------
!!    calling getallq0p returns the following:
      use m_q0p,only: getallq0p,q0i,wt,nq0i,nq0itrue,nq0iadd
!!    calling getbzdata1 returns the following:
      use m_get_bzdata1,only:  getbzdata1,
     &  nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm,
     &  qbz,wbz,qibz,wibz,
     &  qbzw,
     &  idtetf, ib1bz, idteti,
     &  irk, nstar, nstbz,
     &  qbzm, qbzwm
      use keyvalue,only: getkeyvalue
      implicit none
      integer ::nnn(3),ifiqg,ifiqgc,ngcxx,nx0,
     &     ngrp,i,j,iq,iq00,ngp,ngpmx,ngc,ngcmx,nqnum,iq0pin,
     &     nline,nlinemax,ifsyml,iqq,is,nk,ix,nqnumx,i1,ifkpt
      real(8)  :: plat(3,3),qlat(3,3),q(3),qp(3),
     &     QpGcut_psi, QpGcut_Cou,QpGcut,alpv(3),xn,alp,
     &     voltot,pi,alat,tripl,
     &     symops(3,3,ngrp),xx,qqx(3),alpm
      integer,allocatable:: ngvecp(:,:), ngvecc(:,:),
     &     ngpn(:),ngcn(:),ngvect(:,:,:),ngcx(:), nqq(:)
      real(8),allocatable ::
     &     qq(:,:),qq1(:,:),qq2(:,:),qqm(:,:)
      real(8) :: ginv(3,3),aaa
      integer :: ifiqibz,iqibz,mtet(3)
      logical ::tetrai,tetraf,tetra_hsfp0
      integer :: ifbz
!     integer:: bzcase
c     logical :: readgwinput
      integer:: nqnumm,ifiqmtet,verbose,q0pchoice,ifiqbz,iqbz !,auxfunq0p
      real(8)::aaij,bbij
      logical:: qbzreg

      logical :: qreduce ,qreduce0
      real(8),allocatable:: qsave(:,:)
      integer:: imx,ifinin,il,imx0,imxc,imx0c,nnn3(3),imx11(1,1)
      integer,allocatable :: ngvecprev(:,:,:),ngveccrev(:,:,:)

      real(8):: ddq(3)
      logical :: offmesh=.false. ,offmeshg=.false.
      logical :: regmesh=.false. ,regmeshg=.false. ,  timereversal

!     logical ::anyq,newaniso,debug=.true.
      logical ::debug=.true.
!     real(8),allocatable:: qany(:,:)
!     integer:: nany,ifqpnt,ret,imxc,nnn3(3),imx0c,imx11(1,1)
!     real(8):: deltaq,delta5,delta8,deltaq_scale!=1d0/3.0**.5d0

      integer:: nqi,ig,iq0i
      real(8),allocatable:: wti(:),qi(:,:)
      integer:: iclose,iopen
!     real(8):: rrr(3),r2s,emat(3,3),qxx(3)

C     real(8),allocatable:: funa(:,:),wsumau(:),yll(:,:)
C     real(8)::qg(3),alpqg2,qg2,tpiba
      real(8)::volinv,wtrue00,qg(3),alpqg2,qg2,tpiba
      character*99:: q0pf        !nov2012
      integer:: dummyia(1,1),iimx,irradd,nmax
      real(8):: epstol=1d-8,tolq=1d-8,qx(3),qxx(3)
      logical :: newoffsetG !july2014
      real(8),allocatable:: wt0(:)
      integer,allocatable::irr(:)
      real(8):: dq_(3),qlatbz(3,3)
      integer:: gammacellctrl,nnng(3),ifile_handle,ifi0,itet
      real(8)::imat33(3,3)
      logical:: lnq0iadd

c------------------------------------------------
      print *,' mkqg2: '
      qreduce0 = qreduce()
      newoffsetG=.true. !newaniso()
      if(iq0pin == 101) then
         iq0pin=1
         newoffsetG=.false. !for old oldset Gamma case
      endif
!!     band case --- iq0pin == 3 ==>read syml file. E.g for Imag-part calcualtion along a symmetry line.
!!     nqq(is),qq1(1:3,is),qq2(1:3,is),is =1,nline
      if(iq0pin == 3) then
        qreduce0=.false.
        nlinemax = 50
         allocate(nqq(nlinemax),qq1(1:3,nlinemax),qq2(1:3,nlinemax))
         ifsyml = ifile_handle()
         open(ifsyml,file='SYML')
         nline = 0
         do
            nline = nline + 1
            read(ifsyml,*,err=601,end=601)
     &        nqq(nline),qq1(1:3,nline),qq2(1:3,nline)
         enddo
 601     continue
         close(ifsyml)
         nline = nline - 1
         write(6,"(/' Symmetry lines:'/' points',12x,'start',22x,'end')")
         do is=1,nline
            write(6,"(i6,2x,3f8.4,2x,3f8.4)")
     &       nqq(is),(qq1(i,is),i=1,3),(qq2(i,is),i=1,3)
         enddo
         nqnumx = sum(nqq(1:nline))
         allocate( qq(1:3,nqnumx),irr(nqnumx) )
         iqq = 0
         do is = 1,nline
            nk = nqq(is)
            do iq=1,nk
               xx = 0d0
               if(nk>1) xx=(iq-1d0)/(nk-1d0)
               qqx = xx*qq2(1:3,is)+(1d0-xx)*qq1(1:3,is)
               iqq = iqq + 1
               qq(1:3,iqq) = qqx
               write (6,"('  q=',3f7.3)") qq(1:3,iqq)
            enddo
         enddo
         nqnum = iqq
         write (6,"(' Total number of q-points:',i5/)") nqnum
         call minv33tp(plat,qlat) !it was dinv33(plat,1,qlat) by Ferdi
         goto 2001
      endif

!! --- Ordinary case --- iq0pin == 1 or 2
      call getkeyvalue("GWinput","delta",aaa)
      if(aaa<0d0) then
         print * ,'GWinput delta<0: tetrahedron method for x0'
         tetraf=.true.
      else
         print * ,'GWinput delta>0: not use tetrahedron method for x0'
         tetraf=.false.
      endif

      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
C      call dinv33x (plat,qlat)  !it was dinv33(plat,1,qlat) by Ferdi
      call minv33tp(plat,qlat)
      call minv33(qlat,ginv)
      imat33=0d0
      imat33(1,1)=1d0
      imat33(2,2)=1d0
      imat33(3,3)=1d0
      if(sum(abs(matmul(transpose(qlat),plat)-imat33))>tolq) call rx('qlat*plat err')
      if(sum(abs(matmul(ginv,qlat)-imat33))>tolq)            call rx('ginv=qlat^-1 err')
      write(6,*)'=== plat ==='
      write(6,"(3d23.15)") plat
      write(6,*)'=== qlat ==='
      write(6,"(3d23.15)") qlat
C      write(6,*)'=== ginv==='
C      write(6,"(3f9.4)") ginv

      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
      if (sum(abs(mtet))<3) then
         print *, ' we use regular meshing for tetrahedron scheme '
      endif

!! For gammacellctrl==2, we only consider tetrahedron method within the Gammacell.
!! The Gammacell is a part of BZ made from three vectors following qlatbz=(qlat(:,1)/n1q,...)
!! Then the Gamma point is in the middle of micro_qlat = (qlat(:,1)/n1q,qlat(:,2)/n2q,...)
!! To get qbz which is in the Gamma cell, we use shift in the getbzdata1 for gammacellctrl=2.
!! Tetrahedron method is applied for such qbz.
      if(gammacellctrl==2) then
        do i=1,3
          qlatbz(:,i) = qlat(:,i)/nnn(i) !qlat for Gamma cell
        enddo
        call getkeyvalue("GWinput","GammaDivn1n2n3",nnng,3)
        nnn = nnng          !division of Gamma cell
        dq_ = -matmul(qlatbz(1:3,1:3),(/.5d0,.5d0,.5d0/))
!This shift vector is to make the Gamma point centered in the Gamma cell.
        tetrai=.false.
        call minv33(qlatbz,ginv)
        write(6,*)'=== Gammacell qlatgz ==='
        write(6,"(3d23.15)") qlatbz
        write(6,*)'=== Gammacell ginv ==='
        write(6,"(3f9.4)") ginv
cccccccccccccccccccccc
c        qreduce0=.false.
cccccccccccccccccccccc
      else
        qlatbz(:,:) = qlat(:,:)
        tetrai = .true.         !used in heftet tetra_hsfp0()
        dq_ = 0d0
        if(.not.qbzreg()) dq_ = -matmul(qlat(1:3,1:3),(/.5d0/nnn(1),.5d0/nnn(2),.5d0/nnn(3)/))
                          !This dq_ is off-gamma mesh, used when qbzreg=F
      endif
ccccccccccccccccccccccccccc
c      dq_=0d0
ccccccccccccccccccccccccccc
      if(sum(abs(dq_))>tolq) write(6,'(" Shift vector (skip Gamma) by dq_=",3f9.4)')dq_

!! Get BZ data by 'call getbzdata1'
!! See following data after 'use getbzdata1' at the top of this routine.
!! In the case of gammacellctrl=2, we only calculate quantities in the Gamma cell.
!! Thus we have special meanings of nqbz. GWinput --> GammaDivn1n2n3 4 4 4
      call getbzdata1(qlatbz,nnn, !plat bzcase,
     & symops,ngrp,tetrai,tetraf,mtet,gammacellctrl) !all are inputs. output: See use.

      print *,' Writing BZDATA...'
      ifbz = ifile_handle()
      open (ifbz, file='BZDATA')
      write(ifbz,"(10i10)")  nqbz,nqibz, nqbzw, ntetf, nteti,ngrp
      write(ifbz,"(10i10)")  nnn(1:3) !n1q,n2q,n3q

      call rwbzdata(ifbz,-1,ngrp,qlatbz,ginv,
     i    qbz,wbz,nstbz,nqbz,
     i    qibz,wibz,nstar,irk,nqibz,
     i    idtetf,ntetf,qbzw,ib1bz,nqbzw,
     i    idteti,nteti,dq_)
      close(ifbz)

C ... Write QIBZ
      write(6,*)' qibz are written in QIBZ file...'
      ifiqibz = ifile_handle()
      open (ifiqibz,file='QIBZ') !write q-points in IBZ.
      write(ifiqibz,"(i10)") nqibz
      do iqibz = 1,nqibz
         write(ifiqibz,"(3d24.16,3x,d24.16)") qibz(1:3,iqibz),wibz(iqibz)
      enddo
      close(ifiqibz)
C ... Write QBZ
      ifiqbz = ifile_handle()
      open (ifiqbz,file='QBZ') !write q-points in IBZ.
      write(ifiqbz,"(i10)") nqbz
      do iqbz = 1,nqbz
         write(ifiqbz,"(3d24.16,3x,d24.16)") qbz(1:3,iqbz)
      enddo
      close(ifiqbz)
!!  Write KPNTin1BZ.mkqg.chk (files *.chk is only for check.).
      ifkpt = ifile_handle()
      open(ifkpt,file='KPTin1BZ.mkqg.chk')
      write(ifkpt,*)"  qbz --> shorten(qbz)"
      do      i1 = 1,nqbz
         call shorbz(qbz(1,i1),qp,qlat,plat)
         write (ifkpt,"(1x,i7,4f10.5,'   ',3f10.5)")
     &        i1,qbz(1,i1),qbz(2,i1),qbz(3,i1),wbz(i1),qp
      end do
      close (ifkpt)
      write(6,*) ' --- TOTAL num of q =',nqbz
      write(6,*)
      write(6,"( '  ngrp  = ',i3)")ngrp
      write(6,'("  qibz=",i6,3f12.5)')(i,qibz(1:3,i),i=1,min(10,nqibz))
      write(6,*)" ... QIBZ is written in QIBZ file ..."

!!   alpha is for auxiliary function for offset Gamma method.
      call getkeyvalue("GWinput","alpha_OffG",alp,default=-1d60)
      alpv(:)=alp
      if(alp==-1d60) then
        call getkeyvalue("GWinput","alpha_OffG_vec",alpv,3,default=(/-1d50,0d0,0d0/))
        if(alpv(1)==-1d50) then
          call rx( ' mkqg: No alpha_offG nor alpha_offG_vec given in GWinput')
        endif
      endif
      print *
      print *,' alpv=',alpv
      print *
      alpm = minval(alpv)
      if(alpm<=0d0) call rx( 'alpha_offG or alpha_offG_vec <=0')

!! Determine G vectors for q points set by getgv2
      if(iq0pin==1) then        ! get q0x (offset q=0 point)
!! I now think this QpGcut is large enough.
        QpGcut = sqrt(25d0/alpm) !a.u. !exp( -alp*QpGcut**2) !alp * QpGcut**2 = 22
        allocate( ngcx(nqbz) )
        ngcx=1
        do iq = 1, nqbz
          q = qbz(1:3,iq)
          call getgv2(alat,plat,qlat,q, QpGcut, 1, ngcx(iq),  dummyia)
        enddo
        ngcxx = maxval(ngcx)
        allocate( ngvect(3,ngcxx,nqbz) )
        do iq = 1, nqbz
          q = qbz(1:3,iq)
          call getgv2(alat,plat,qlat, q, QpGcut, 2,
     &      ngcx(iq), ngvect(1:3,1:ngcx(iq),iq) )
        enddo
      endif


!! getallq0p all inputs
!! Q0P is offset Gamma or k point given in GWinput
!! see use m_q0p =>  q0i,wt,nq0i,nq0itrue are outputs
!! we now have q0i(:,nq0i+1,nq0i+nq0iadd).
!!    q0i(:,1:nq0i+n0qiadd) contains all q0x(:,i)= qlat(:,i)/nnn(i)/2d0*deltaq_scale() for i=1,3.
c      lnq0iadd=.true.
      call getallq0p(iq0pin,newoffsetG,alat,plat,qlat,nnn,alp,alpv, !apr2016
     i ngcxx,ngcx,nqbz,nqibz,nstbz,qbz,qibz,symops,ngrp,ngvect,lnq0iadd)
c      print *,'size q0i=',size(q0i),ubound(q0i),lbound(q0i)
      do i=nq0i+1,nq0i+nq0iadd
        write(6,"('  q0iadd=  ', i3, 3f10.5)") i,q0i(:,i)
      enddo

!! Four kinds of mesh points. Q0P means offset Gamma (slightly different from Gamma).
!! Which we need?
!! 1. regular
!! 2. offregular (not including Gamma)
!! 3. regular    + Q0P
!! 4. offregular + Q0P
      if(iq0pin==2) then        !this is just for dielectric case
         regmesh = qbzreg()
      else
         regmesh = .true.
      endif
      regmeshg = qbzreg()       !Gamma mesh based on regular mesh
      offmesh =  .not.qbzreg()  !we fix bzcase=1 now. apr2015.
      offmeshg = .not.qbzreg()  !Gamma mesh based on off-regular mesh
      print *,' regmesh offmeshg=', regmesh,regmeshg !regular,     regular+shifted
      print *,' offmesh offmeshg=', offmesh,offmeshg !offregmesh, offregular+shifted

!!  We check wether all q0i \in qbz or not. <--- Takao think this block is not necessary now.
      call minv33(qlat,ginv)
      nqnum = nqbz
      allocate( qq(1:3,nqnum),irr(nqnum) )
      qq(1:3,1:nqbz) = qbz(1:3,1:nqbz)
      do iq0i=1,nq0i+nq0iadd
         do iq=1,nqbz
            if(sum(abs(q0i(:,iq0i)-qq(:,iq)))<tolq) goto 2112
            call rangedq( matmul(ginv,q0i(:,iq0i)-qq(:,iq)), qx)
            if(sum(abs(qx))< tolq) goto 2112
         enddo
         goto 2111
 2112    continue
         qq(:,iq) = q0i(:,iq0i) !replaced with equivalent q0i.
      enddo
      print *,' --- We find all q0i in qbz. Skip qreduce.'
      goto 2001
 2111 continue

!! Accumulate all required q points
      deallocate(qq,irr)
      nqnum = nqbz  + nqbz*(nq0i+nq0iadd)
      nqnum = nqnum + 1         !add Gamma
      nqnum = nqnum + nq0i + nq0iadd      !add Gamma + q0i
      allocate( qq(1:3,nqnum),irr(nqnum) )
      ix = 0
      if(regmesh) then
         qq(1:3,1:nqbz) = qbz(1:3,1:nqbz)
         ix = ix+ nqbz
      endif
!!  - Off Regular mesh.
      if(offmesh) then
         do iq = 1, nqbz
            ix = ix+1
            qq(1:3,ix) = qbz(1:3,iq) - dq_
         enddo
      endif
c      nnn   =  ix               !n1q*n2q*n3q!      if(offmesh) nnn = 2*n1q*n2q*n3q
c      print *,' nnn=',nnn       !This is the number to calcualte Vxc
!!  - Shifted mesh
      if(regmeshg) then
         do iq00 = 1, nq0i+ nq0iadd
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbz(1:3,iq) +  q0i(1:3,iq00)
            enddo
         enddo
      endif
      if(offmeshg) then
         do iq00 = 1, nq0i+ nq0iadd
            do iq   = 1, nqbz
               ix = ix+1
               qq(1:3,ix) = qbz(1:3,iq) - dq_ + q0i(1:3,iq00)
            enddo
         enddo
      endif
!!  - Add offset Gamma and Gamma point (these can be removed by qreduce and q0irre)
      do iq00 = 1, nq0i+ nq0iadd
         ix = ix+1
         qq(1:3,ix) = q0i(1:3,iq00)
      enddo
      ix=ix+1
      qq(1:3,ix)=0d0


!! (this mtet block is not used now) Get qqm; q point for eigenvalues.
!! Saved to Qmtet. Not so much used now...
!! We need check again if we like to use this branch again (2016apr)
      if(sum(abs(mtet))/=3) then
         nqnumm= nqbzm * (nq0i+ nq0iadd +1)
         allocate( qqm(1:3,nqnumm) )
         ix=0
         do iq00 = 1, 1 + nq0i+ nq0iadd
            do iq   = 1, nqbzm
               ix = ix+1
               if(iq00==1) then
                  qqm(1:3,ix) = qbzm(1:3,iq)
               else
                  qqm(1:3,ix) = q0i(1:3,iq00-1) + qbzm(1:3,iq)
               endif
            enddo
         enddo
         ifiqmtet=ifile_handle()
         open(ifiqmtet, file='Qmtet')
         write(ifiqmtet,"(i10)") nqnumm
         do iq=1,nqnumm
            write(ifiqmtet,"(3d24.16)") qqm(1:3,iq)
         enddo
         close(ifiqmtet)
         deallocate(qqm)
      endif

!! Remove equivalent q point by the translational symmetry
      if( qreduce0 ) then
         print *,'call qqsave with nq0i nq0iadd nqnum',nq0i,nq0iadd,nqnum
         call cputid (0)
         nmax= nq0i+nq0iadd+nqnum
         allocate(qsave(3,nmax)) !,qsavel(nmax))
         imx=0
         if(iq0pin /=1) then
            do iq=1,nq0i+ nq0iadd
               call qqsave(q0i(1:3,iq),nmax,ginv,qsave,imx)
            enddo
         endif
         do iq=1,nqnum
            call qqsave(qq(1:3,iq),nmax,ginv,qsave,imx)
         enddo
         nqnum = imx
         qq(:,1:imx)=qsave(:,1:imx)
         deallocate(qsave)
      endif
!! ------------------------------------------
 2001 continue
!! ------------------------------------------


!! Here we get all required q points. We do reduce them by space group symmetry.
      if(allocated(wt0)) deallocate(wt0)
      allocate(wt0(nqnum+nq0i+ nq0iadd ),qi(3,nqnum+nq0i+ nq0iadd ),wti(nqnum+nq0i+ nq0iadd ))
      wt0=1d0
!! Set irreducible k-point flag. irr=1 for (irredusible point) flag, otherwise =0.
!! irr(iq)=1 for irreducile qq(:,iq), iq=1,nqnum
      call q0irre(qibz,nqibz,qq,wt0,nqnum,symops,ngrp, qi,nqi,wti,plat,.true.,0,irr)
!! nqnum is the finally obtained number of q points.
      allocate(ngpn(nqnum), ngcn(nqnum))
      if(debug) write(6,*) ' --- q vector in 1st BZ + Q0P shift. ngp ---'
      imx=0
      imxc=0
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q
         if(iq0pin==1) then !use qxx on regular mesh points if q is on regular+Q0P(true).
            do iqbz=1,nqbz
            do i=1,nq0itrue+ nq0iadd  ! nq0itrue/=nq0i for anyq=F nov2015
               if(sum(abs(qbz(1:3,iqbz)-dq_+ q0i(:,i)-qxx))<tolq) then
                  qxx=qbz(1:3,iqbz)
                  exit
               endif
            enddo
            enddo
         endif
         ngpn(iq)=1
!! get nqpn. # of G vector for |q+G| < QpGcut_psi
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi,1,ngpn(iq),imx11) !imx11 !nov2015
         imx0=imx11(1,1)
         if(imx0>imx) imx=imx0
         ngcn(iq)=1
!! get ngcn. # ofG vector for |q+G| < QpGcut_cou
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou,1,ngcn(iq),imx11) !imx11 to avoid warning.
         imx0c=imx11(1,1)
         if(imx0c>imxc) imxc=imx0c
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngpn(iq) !,ngcn(iq,iq00)
         if(verbose()>150)write(6,'(3f12.5,3x,2i4)') q ,ngcn(iq) !,ngcn(iq,iq00)
      enddo

!! Get G vectors and Write q+G vectors -----------
      ngpmx = maxval(ngpn)
      ngcmx = maxval(ngcn)
      write(ifiqg ) nqnum,ngpmx,QpGcut_psi,nqbz,nqi,imx,nqibz
      write(ifiqgc) nqnum,ngcmx,QpGcut_cou,nqbz,nqi,imxc
!! :nqi:   The number of irreducible points (including irr. of offset points). irr=1.
!! ::       We calcualte eigenfunction and Vxc for these points.
!! :nqnum: total number of q points.
!! :imx:   to allocate ngvecprev as follows.
      print *,' number of irreducible points nqi=',nqi
      print *,' imx nqnum=',imx,nqnum
      write(6,*) ' --- Max number of G for psi =',ngpmx
      write(6,*) ' --- Max number of G for Cou =',ngcmx
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx) )       !inverse mapping table for ngvecp (psi)
      allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc) ) !inverse mapping table for ngvecc (cou)
      ngvecprev=9999
      ngveccrev=9999
      do iq = 1, nqnum
         q = qq(1:3,iq)
         qxx=q
         q0pf=''
         do iqbz=1,nqbz  !use qxx on regular mesh points if q is on regular+Q0P(true).
         do i=1,nq0itrue+ nq0iadd  !nq0itrue/=nq0i for anyq=F nov2015
            if(sum(abs(qbz(1:3,iqbz)-dq_+ q0i(:,i)-qxx))<tolq) then
               if(sum(abs(q0i(:,i)-qxx))<tolq) then
                 q0pf=' <--Q0P  '   ! offset Gamma points
               else
                 q0pf=' <--Q0P+R'   ! offset Gamma points-shifted nov2015
               endif
               if(iq0pin==1) then
                  qxx=qbz(1:3,iqbz)
               endif
               exit
            endif
         enddo
         enddo
         ngp = ngpn(iq)
         ngc = ngcn(iq)
         write(6,"(' iq=',i8,' q=',3f9.5,' ngp ngc= ',2i6,' irr.=',i2,a)") !irr=1 is irreducible k points.
c         write(6,"(' iq=',i8,' q=',3f17.13,' ngp ngc= ',2i6,' irr.=',i2,a)") !irr=1 is irreducible k points.
     &        iq, q, ngp, ngc, irr(iq),trim(q0pf)
         allocate( ngvecp(3,max(ngp,1)), ngvecc(3,max(ngc,1)) )
         call getgv2(alat,plat,qlat, qxx, QpGcut_psi, 2, ngp,  ngvecp) ! for eigenfunctions (psi)
         call getgv2(alat,plat,qlat, qxx, QpGcut_Cou, 2, ngc,  ngvecc) ! for Coulomb        (cou)
         write (ifiqg) q, ngp, irr(iq)
         do ig = 1,ngp
            nnn3 = ngvecp(1:3, ig)
            ngvecprev( nnn3(1), nnn3(2),nnn3(3)) = ig
         enddo
         write (ifiqg)  ngvecp,ngvecprev !ngvecprev is added on mar2012takao
         do ig = 1,ngc
            nnn3 = ngvecc(1:3, ig)
            ngveccrev( nnn3(1), nnn3(2),nnn3(3)) = ig
         enddo
         write (ifiqgc) q, ngc
         write (ifiqgc) ngvecc,ngveccrev
         deallocate(ngvecp,ngvecc)
      enddo
      deallocate(ngpn,ngcn,ngvecprev,ngveccrev)
      if(iq0pin==1) deallocate(ngvect)
      if(debug) print *,'--- end of mkqg2 ---'
      end

!! #############################################################################
C      doubleprecision function tripl(a,b,c)
C!     ! == tripl  ==
C      implicit real*8 (a-h,p-z), integer(o)
Cc     dimension a(3),b(3),c(3)
Cc     ... g77 needs this rewriting for opt > O1
C      dimension a(0:2),b(0:2),c(0:2)
Cc     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
Cc     .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
Cc     ... g77 needs this rewriting
C      tmp = 0.d0
C      do i = 0,2
C         j = mod(i + 1, 3)
C         k = mod(i + 2, 3)
C         tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
C      enddo
C      tripl = tmp
C      end


!! #############################################################################
C      subroutine genbzdata1_xxxxxxx(plat,n1q,n2q,n3q,symops,ngrp,
C     &     tetraf, tetrai,
C     o     qbz, wbz,              nqbz,
C     o     qibz,wibz, nstar,irk,  nqibz,
C     o     idtetf, ntetf, qbzw,ib1bz, nqbzw,
C     o     idteti, nteti )
C!     ! == genbzdata1_xxxxxxx: Generate BZ mesh data reuired for GW ==
CC     iput
CC     plat,qlat,n1q,n2,n3q,symops,ngrp,
CC     tetraf: if true, this generate the BZ data for full-1stBZ tetrahedron method
CC     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
CC     tetrai: if true, this generate the BZ data for tetrahedron method in irreducible BZ.
CC     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
CC     output
CC     o  qbz, wbz,              nqbz,
CC     o  qibz,wibz, nstar,irk,  nqibz,
CC     ouptput(optional)
CC     o  idtetf, ntetf, qbzw,ib1bz, nqbzw
CC     o  idteti, nteti
CC-------------------------------------------------------------------------
C      implicit none
C      integer:: nqbz,n1q,n2q,n3q,ntetf,nteti,nqbzw,iqbz,ifbz
C     &     ,nqibz,iqibz,itet,ngrp
C      real(8) :: plat(3,3),qbas(3,3),ginv(3,3),vol,symops(3,3,ngrp)
C     &     ,qbasmc(3,3)
C      logical tetraf,tetrai
C      integer,allocatable:: ipq(:),iw1(:)
C      real(8):: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:)
C      integer :: idtetf(:,:),ib1bz(:),idteti(:,:),irk(:,:)
C     &     ,nstar(:)
C      end
C
C
C!! #############################################################################
      integer function llxxx(ilm)
      integer,parameter :: lmx=50
      integer,save:: lla((lmx+1)**2)
      logical:: init=.true.
      if (ilm>(lmx+1)**2) call rx( 'll: ilm too large')
      if (init) then
        do l=0,lmx
          lini= l**2 + 1
          lend=(l+1)**2
          lla(lini:lend)=l
        enddo
      endif
      llxxx = lla(ilm)
      return
      end
