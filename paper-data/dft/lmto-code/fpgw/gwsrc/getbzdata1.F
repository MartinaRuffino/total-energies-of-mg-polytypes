      real(8) function xqconv(x)
!!! x --> xqconv : uniform to non uniform mesh converter
!! x is [0,1] --> xqcon = [0,1]
!! x can be -1 <= x =< 1
      use keyvalue,only: getkeyvalue
      real(8),intent(in)::x
      real(8),parameter:: pi=3.1415926535897932d0
      real(8),save:: adiv,bdiv
      logical,save:: oncew=.true.
      logical:: ggg
!! BZ division setting.
      if(oncew) then
         inquire(file='GWinput',exist=ggg)
         if(.not.ggg) then
            adiv=1d0
         else
            call getkeyvalue("GWinput","BZadiv",adiv,default=1d0)
         endif
         write(6,"('BZadiv= ',f6.3)") adiv
         bdiv = (1d0 -adiv)/2d0
         oncew=.false.
      endif
!!   adiv+2b=1
      xqconv           = adiv*x + bdiv*(1-cos(pi*x))
      if(x>1d0) xqconv = adiv*x + bdiv*(1-cos(pi*(x-1d0))+2d0)
      if(x<0d0) xqconv = adiv*x + bdiv*(1-cos(pi*(x+1d0))-2d0)
c      print *,'x xqconv=',x,xqconv
      end function xqconv

      subroutine xconvv(xin,xout)
!!! xin is converted to xout (uniformmesh to non-uniform mesh).
      real(8)::xin(3),xout(3),xqconv
      integer:: i
      do i=1,3
        xout(i)=xqconv(xin(i))
      enddo
      end subroutine xconvv

      module m_get_bzdata1
C- Make k-points
C  In addition to module variables below, this also write mtet file (search ifmtet) for mtet mode.
c------------------------------------------------------------
      implicit none  ! all are outputs
      integer,protected::             nqbz, nqibz, nqbzw,ntetf,nteti,nqbzm
      real(8),allocatable,protected:: qbz(:,:),wbz(:),qibz(:,:),wibz(:),qbzw(:,:)
      real(8),protected::             dq_(3)!,qbasmc(3,3) !,dq_bzcase2(3),
      integer,allocatable,protected:: idtetf(:,:),ib1bz(:),idteti(:,:),irk(:,:),nstar(:),nstbz(:)
      real(8),allocatable,protected:: qbzm(:,:),qbzwm(:,:)

c------------------------------------------------------------
      contains
      subroutine getbzdata1(qlat,nnn,
     &    symops,ngrp,tetrai,tetraf,mtet,gammacellctrl)
      use keyvalue,only: getkeyvalue
      implicit none
      logical,intent(in) :: tetrai,tetraf
      integer,intent(in) :: ngrp,nnn(3),mtet(3),gammacellctrl
      real(8),intent(in) :: qlat(3,3),symops(3,3,ngrp)
C Local
      integer,allocatable:: IPQ(:,:,:),indexkw(:,:,:),indexk(:,:,:)
      real(8),allocatable:: qcm(:,:,:),wtet(:,:,:)
      logical:: skipgammacell
      integer,allocatable :: idtetfm(:,:,:),ib1bzm(:) !multipled tetrahedron.
      real(8):: qbzshift(3)
      integer:: nmtet,nqbzwm,ntetfm,ifmtet,n1qm,n2qm,n3qm,nnnx(3),itet,ix,im
      integer:: iccc,nadd=0
      real(8):: plat(3,3),ginv(3,3),qc(3,0:3),qbi(3,3),qmic(3,3)
      integer:: nnnv(3),ngcell
      real(8), parameter:: tolq=1d-6
      procedure(logical) :: qbzreg

c-------------------------------------------------------
      call minv33tp (qlat,plat)  !qlat --> plat
      call minv33(qlat,ginv)

      qbzshift=0d0
c     if(icase==2) qbzshift=0.5d0
c     if(icase>2) call rx( 'mkqg: wrong bzcase')

      if(gammacellctrl==2) then
        nadd=1  ! nddd=1 give end point of BZ (for gammacell case). See nqbz
        qbzshift=-1d0/2d0       ! gamma centered.
              !this was (/n1q/2d0,n2q/2d0,n3q/2d0/) !this shift is so that gamma point is centerd.
      endif

      nnnv = nnn+nadd
c      wfac=dble(product(nnn))/product(nnnv)
      nqbz = product(nnnv)      !(nnn(1)+nadd)*(nnn(2)+nadd)*(nnn(3)+nadd)
      allocate(qbz(3,nqbz),wbz(nqbz)) !,nstbz(nqbz))
      allocate(nstbz(nqbz)
     & ,indexk (0:nnnv(1)-1,0:nnnv(2)-1,0:nnnv(3)-1)
     & ,indexkw(0:nnn(1),0:nnn(2),0:nnn(3)))
      nstbz(1:nqbz)=0

      write(6,"(' getbzdata1: n1n2n3 =      ',3i5)") nnn
      write(6,"(' getbzdata1: n1n2n3 +nadd =',3i5)") nnnv


      call genqbz(qlat,nnn(1),nnn(2),nnn(3),qbz,wbz,nstbz,nadd,qbzshift)
C      allocate(qbzz(3,nnnv(1),nnnv(2),nnnv(3)),qbzw(nnnv(1)+1,nnnv(2)+1,nnnv(3)+1)))
C      weight = 1d0/dble(nnn1t*nnn2t*nnn3t)
C
      allocate(qibz(3,nqbz),ipq(nnnv(1),nnnv(2),nnnv(3)),wibz(nqbz))
      call bzmesh(plat,qmic,nnn(1),nnn(2),nnn(3),symops,ngrp,ipq,qibz,wibz,nqibz,nqbz,nadd,qbzshift) !Make q-points in IBZ.

      allocate(nstar(nqibz),irk(nqibz,ngrp))
      call nkstar  (qibz,qbz,symops,ginv,
     d              nqibz,nqbz,ngrp,gammacellctrl,
     o              nstar,irk )

      ntetf=-1
      nteti=-1
      nqbzw=-1
      if (tetraf) then
        nqbzw = (nnn(1)+1)*(nnn(2)+1)*(nnn(3)+1)
        ntetf = 6*nqbz
        ngcell = 1
        allocate( idtetf(0:3,ntetf), ib1bz(nqbzw), qbzw(3,nqbzw) )
        skipgammacell = .false.
        if(gammacellctrl==1) skipgammacell=.true.
C       call tetfbzf(qlat,nnn(1),nnn(2),nnn(3),qbz,nqbz,skipgammacell,1,qbzshift,nadd,idtetf,qbzw,ib1bz,ntetf)
        call tetfbz(qlat,nnn(1),nnn(2),nnn(3),nadd,skipgammacell,ngcell,qbzshift,1d0,idtetf,qbzw,ib1bz,iccc)
        if (iccc /= ntetf) call rx('problem with tetfbz')
        print *, 'debugging sum idtetf, qbzw,ib1bz',sum(idtetf),sum(qbzw),sum(ib1bz)
      endif
      if(tetrai) then
        allocate(idteti(0:4,6*nqbz))
        qmic(:,1)= qlat(:,1)/dble(nnn(1))
        qmic(:,2)= qlat(:,2)/dble(nnn(2))
        qmic(:,3)= qlat(:,3)/dble(nnn(3))
        call tetirr(qmic,nnn(1),nnn(2),nnn(3),ipq,nqibz,nteti,idteti)
        print *, 'debugging sum idteti',sum(idteti)
      endif
      deallocate(ipq)

c     ---for multipled tetrahedron-----------------
      if(sum(abs(mtet))/=3) then
        print *, 'multitet mode: mtet=',mtet
        n1qm = mtet(1)*nnn(1)
        n2qm = mtet(2)*nnn(2)
        n3qm = mtet(3)*nnn(3)
        nmtet   = mtet(1)*mtet(2)*mtet(3)
        nqbzm   = nmtet * nqbz
        nqbzwm  = (n1qm+1)* (n2qm+1)* (n3qm+1)
        ntetfm  = ntetf * nmtet
        allocate(
     &       idtetfm(0:3,nmtet,ntetf), qbzwm(3,nqbzwm),
        ! Index for tetrahedron;    qbzmw(idtetfm) gives extended q vector.
     &       ib1bzm(nqbzwm), qbzm(3,nqbzm) )
        ! qbzm(1:3,ib1bz(iq)) gives q vector within the 1st bz.
! The datas idetetfm, qbzmw, ib1bzm, nmete, ntetf, nqbzm,
!     eigen(nband,nqbzm) are required for the multiply-divided tetrahedron method.
!
        allocate( qcm(1:3,0:3, nmtet), wtet(0:3, nmtet, ntetf),
     &          indexkw(0:n1qm,0:n2qm,0:n3qm) )
        qmic(:,1)= qlat(:,1)/dble(nnn(1)*mtet(1))
        qmic(:,2)= qlat(:,2)/dble(nnn(2)*mtet(2))
        qmic(:,3)= qlat(:,3)/dble(nnn(3)*mtet(3))
        call minv33(qmic,qbi)
!       qbi = transpose(qbi)  ? there, or not?
        call qwider(1,qmic, mtet(1)*nnn(1),mtet(2)*nnn(2),mtet(3)*nnn(3),
     &                       n1qm,n2qm,n3qm, ib1bzm, qbzwm, qbzm, indexkw )
        do itet=1, ntetf  !--- idtetfm(0:3, nmtet, ntetf), nmtet divition.
          do ix=0,3
            qc(:,ix)= qbzw(:, idtetf(ix,itet))
          enddo
ccccccccccccccccccccccccccccccccccccccccccccccc
c        do ix = 1,3
c        kvec(1:3,ix) = qc(1:3,ix) - qc(1:3,0)
c        enddo
c        write(6,"('itet vol=',i5,d13.5)") itet,abs(det33(kvec(1:3,1:3))/6d0)
cccccccccccccccccccccccccccccccccccccccccccccccc
          call tetdevide(qc,qcm, wtet(:,:,itet),mtet(1),mtet(2),mtet(3))   !qc ---> qcm
          do im=1,nmtet ! Index for micro-devided tetrahedron
cccccccccccccccccccccccccccccccccccccccccccccccc
c        do ix = 1,3
c        kvec(1:3,ix) = qcm(1:3,ix,im) - qcm(1:3,0,im)
c        enddo
c        write(6,"('itet im vol=',2i5,d13.5)") itet,im,abs(det33(kvec(1:3,1:3))/6d0)
cccccccccccccccccccccccccccccccccccccccccccccccc
            do ix=0,3     ! its four corners
              !qcm ---> idtetfm(ix,im,itet)
              nnnx(1:3) = matmul (qbi,qcm(1:3,ix,im))+ 1d-10
              idtetfm(ix,im,itet) = indexkw(nnnx(1),nnnx(2),nnnx(3))
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                write(6,"(' nnnx=',3i3)") nnnx(1:3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ! qbzm is given as qbzm(1:3, ib1bz(idtetfm(ix,im,itet)) )
              if(abs(sum(nnnx - matmul (qbi,qcm(1:3,ix,im))))>tolq) then
                call rx( 'getbzdata1: nnn is not integer')
              endif
c               write(6,"('n1 n2 n3=',3i5)") nnnx(1:3)
            enddo
          enddo
        enddo
c--- write mtet
        ifmtet=501
        open (ifmtet, file='mtet',form='unformatted')
        write(ifmtet) nmtet,nqbzwm,nqbzm,ntetfm !,n_index_qbzm
        write(ifmtet) idtetfm,ib1bzm,qbzm,qbzwm,wtet  !,index_qbzm
        close(ifmtet)
      endif

      end subroutine getbzdata1

      subroutine nkstar (qibz,qbz,grp,ginv,
     d nqibz,nqbz,ngrp,gammacellctrl,
     o nstar,irotk)
c 92.02.22
c generates the no. stars of k
c i.e. the no. times k appears in the FBZ
c qibz  = k { IBZ
c qbz   = k { FBZ
c grp   = rotation matrices
c nqibz = no. k { IBZ
c nqbz  = no. k { FBZ
c ngrp  = no. rotation matrices
c nstar(k) = no. times k appears in the FBZ
c irotk(k{IBZ,R) = index to k{FBZ
      implicit none
      integer:: nqibz,nqbz,ngrp,ir,k,kp,nsum,ivsum
      real(8):: qibz(3,nqibz),qbz(3,nqbz),grp(3,3,ngrp),ginv(3,3)
      integer:: nstar(nqibz),irotk(nqibz,ngrp)
      integer:: verbose,kout,nirotk,gammacellctrl
      real(8):: diff(3),diff2(3),tolq=1d-6
      if(verbose()>104) then
        print *,' nkstar:'
        do kp = 1,nqbz
          write(6,"(' === kp=',i8,' qbz=',3f8.3)")kp,qbz(:,kp)
        enddo
      endif
      irotk=0
      nstar=0
      do kp = 1,nqbz
        do k  = 1,nqibz
          do ir = 1,ngrp
            if(verbose()>104) print *,' grp=',ir !;      print *, grp(:,ir)
            diff = matmul(grp(:,:,ir),qibz(:,k)) - qbz(:,kp)
            if(gammacellctrl/=2) then
              call rangedq(matmul(ginv,diff), diff2)
            else
              diff2=diff
            endif
            if(verbose()>104) write(6,"(' matmul(ginv,diff)=',3f8.3,' ',3f8.3)") diff, matmul(ginv,diff)
            if(sum(abs(diff2))< tolq) then
              irotk(k,ir)= kp
              kout=k
              nstar(k)   = nstar(k) + 1
              goto 1022
            endif
          enddo
        enddo
        call rx( 'nkstar: can not find irotk')
 1022   continue
c        write(6,"('   kp=',i8,' qbz=',3f18.14, 'ibz qibz',i8,3f18.14)")kp,qbz(:,kp),kout,qibz(:,kout)
      enddo
c      do k  = 1,nqibz
c          write(6,"(' k nstar=',3i6)") k,nstar(k)
c      enddo
c      print *,'sum nstart=',sum(nstar)
      nirotk=0
      do k=1,nqibz
      do ir=1,ngrp
         if(irotk(k,ir)/=0) then
            nirotk=nirotk+1
         endif
      enddo
      enddo
ccccccccccc

!! check that the sum of stars equal to the no. k{FBZ
      nsum = ivsum (nstar,nqibz)
      if (nsum .ne. nqbz) then
        print *,' nums nqbz=',nsum,nqbz
        do k  = 1,nqibz
          do ir = 1,ngrp
            if(irotk(k,ir)/=0) write(6,"(' k ir irotk=',3i6)") k,ir, irotk(k,ir)
          enddo
        enddo
        call rx( 'nkstar: wrong no. stars')
      endif
c$$$c write k { IBZ and no. stars to file KPNT
c$$$      ifkpnt     = ifile('KPNT')
c$$$      if (ifkpnt .gt. 0) then
c$$$        write (ifkpnt,*) 'irreducible k-points and no. stars'
c$$$        write (ifkpnt,*) 'k, k-vector, nstar '
c$$$        do       k = 1,nqibz
c$$$          write(ifkpnt,"(1x,i5,3f8.5,i3)")k,qibz(1,k),qibz(2,k),qibz(3,k),nstar(k)
c$$$        end do
c$$$      endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c      stop 'test end'
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print *,'number irotk and sum(irotk)',nirotk,sum(irotk)
      end subroutine nkstar
      end module

c--------------------------------------------------------------------------
      subroutine qwider(icase,qb,n1,n2,n3,n1w,n2w,n3w, ib1bz,qbzw,qbz,indexkw)
!! == Wider q point mesh. ==
      implicit none
      integer(4):: i1,i2,i3, kount,icase,n1w,n2w,n3w
      integer(4):: n1,n2,n3,ib1bz(*),indexkw(0:n1w,0:n2w,0:n3w)
      integer(4):: indexk(0:n1-1,0:n2-1,0:n3-1)
      real(8):: qb(3,3),qbzw(1:3,*),qbz(1:3,*),hf
      hf=0d0
c      if(icase==2) hf=0.5d0
      kount      = 0
      do      i1 = 1,n1
        do      i2 = 1,n2
          do      i3 = 1,n3
            kount    = kount + 1
            indexk(i1-1,i2-1,i3-1) = kount
            qbz(1:3,kount) = qb(1:3,1)*(i1-1+hf) +qb(1:3,2)*(i2-1+hf) +qb(1:3,3)*(i3-1+hf)
          end do
        end do
      end do
c
      kount      = 0
      do      i1 = 1,n1w+1
        do      i2 = 1,n2w+1
          do      i3 = 1,n3w+1
            kount    = kount + 1
            indexkw(i1-1,i2-1,i3-1) = kount
            qbzw(1:3,kount) =
     &    qb(1:3,1)*(i1-1+hf) + qb(1:3,2)*(i2-1+hf) + qb(1:3,3)*(i3-1+hf)
            ib1bz(kount) = indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
          end do
        end do
      end do
      end
!!--------------------
      subroutine tetdevide(qc, qcm, wtet, mt1,mt2,mt3)   !qc ---> qcm
      integer(4):: mt1,mt2,mt3
      real(8):: qc(3,0:3), qcm(3,0:3,mt1*mt2*mt3),qq(3,0:9),wt(0:3,0:9)
     &         ,wtet(0:3,mt1*mt2*mt3)
      integer(4):: iq(0:3,8),itet,ic
!! == four tetrahedrons at ends ==
      iq(:,1) = (/0,7,8,9/)
      iq(:,2) = (/1,5,6,7/)
      iq(:,3) = (/2,4,6,9/)
      iq(:,4) = (/3,4,5,8/)
c octahedron into four tetrahedron.
      iq(:,5) = (/5,9,6,7/)
      iq(:,6) = (/5,9,7,8/)
      iq(:,7) = (/5,9,8,4/)
      iq(:,8) = (/5,9,4,6/)
      if(mt1==2.and.mt2==2.and.mt3==2) then
        qq(:,0:3)  =  qc(:,0:3)
        qq(:,3+1)  = (qc(:,2) +qc(:,3))/2d0
        qq(:,3+2)  = (qc(:,3) +qc(:,1))/2d0
        qq(:,3+3)  = (qc(:,1) +qc(:,2))/2d0
        qq(:,6+1)  = (qc(:,1) +qc(:,0))/2d0
        qq(:,6+2)  = (qc(:,3) +qc(:,0))/2d0
        qq(:,6+3)  = (qc(:,2) +qc(:,0))/2d0
c
        wt(:,0)  =  (/1d0,0d0,0d0,0d0/)
        wt(:,1)  =  (/0d0,1d0,0d0,0d0/)
        wt(:,2)  =  (/0d0,0d0,1d0,0d0/)
        wt(:,3)  =  (/0d0,0d0,0d0,1d0/)
        wt(:,3+1)  = (wt(:,2) +wt(:,3))/2d0
        wt(:,3+2)  = (wt(:,3) +wt(:,1))/2d0
        wt(:,3+3)  = (wt(:,1) +wt(:,2))/2d0
        wt(:,6+1)  = (wt(:,1) +wt(:,0))/2d0
        wt(:,6+2)  = (wt(:,3) +wt(:,0))/2d0
        wt(:,6+3)  = (wt(:,2) +wt(:,0))/2d0
        do itet=1,8
          do ic=0,3
            qcm (:,ic,itet) = qq(:,iq(ic,itet))
          enddo
          wtet(0,itet) =  sum( wt(0,iq(:,itet)) )/4d0
          wtet(1,itet) =  sum( wt(1,iq(:,itet)) )/4d0
          wtet(2,itet) =  sum( wt(2,iq(:,itet)) )/4d0
          wtet(3,itet) =  sum( wt(3,iq(:,itet)) )/4d0
c         write(6,"(' itet wtet=',i5,5f8.3)")itet,wtet(:,itet),sum(wtet(:,itet))
          if(abs(sum(wtet(:,itet))-1d0)>tolq) call rx( 'tetdevide: sumwtet/=1')
        enddo
      else
        call rx( ' tetdvide: only 2 2 2 has already implimented.')
      endif
      end
