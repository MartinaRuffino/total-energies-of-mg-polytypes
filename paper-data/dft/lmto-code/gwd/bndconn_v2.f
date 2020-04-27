      subroutine bndconn_v2(alat,plat,qlat,
     &      nbas, rmt, bas, ndimh, ldim2,
     &      evl,  ngp, ngp_p,  ngvecp,ngvecp_p,   geig,geig_p,
     &      cphi, cphi_p, ovvx, !ovvx is for non-orthogonal phi !May 2002
     o iiyf,ovv)
c-Get the connenctivity of eigenfunctions along the q points.
co iiyf: connection. idim'th band is mapped to iiyf(idim) band.
co ovv : information of a deformed overlap matrix.
      implicit none
      integer:: ndimh,ldim2,nbas,ngp,ngp_p,iix,iiy,iyiy,istart,j
      integer :: iiyf(ndimh),ifdebug
     &  ,ngvecp(3,ngp),ngvecp_p(3,ngp_p),ifill(ndimh)
      real(8) ::  eigdeg, sumdeg,alat,plat(3,3),qlat(3,3)
     &           ,ovv(ndimh,ndimh),rmt(nbas),bas(3,nbas),evl(ndimh)
      complex(8):: cphi(ldim2,ndimh), cphi_p(ldim2,ndimh)
     & ,geig(ngp, ndimh),geig_p(ngp_p, ndimh)
      complex(8),allocatable::  gpg(:,:),ppovl2(:,:)
      real(8) :: ovvx(ldim2,ldim2)
      complex(8),allocatable:: oc(:,:),coc(:,:),cc(:,:),ovvxc(:,:)
c----
      ifdebug = 1198
c      print *, ' bndconn_v2',cpusec()
c      write(ifdebug,*) ' test xxxxxxxxxxxxxxxxx'
      allocate(gpg(ndimh,ndimh), ppovl2(ngp,ngp_p))
c      print * ,' goto mkppovl2'
c     & ,sum(abs(geig)),sum(abs(geig_p))
c--- Get the connection. ovelap matrix of eigenfunctions
      call mkppovl2(alat,plat,qlat,
     &      ngp,   ngvecp,
     &      ngp_p, ngvecp_p,
     &      nbas, rmt, bas,
     o      ppovl2)
c      print * ,' end of mkppovl2'
c     & ,sum(abs(geig)),sum(abs(geig_p)),sum(abs(ppovl2))
      gpg = matmul(dconjg(transpose(geig)),matmul(ppovl2,geig_p))
      write(ifdebug,'("     ",255i4)') (j,j=1,ndimh)

c---------
      allocate(oc(ldim2,ndimh),coc(ndimh,ndimh),cc(ndimh,ldim2)
     &         ,ovvxc(ldim2,ldim2))
      ovvxc = ovvx
      call matm(ovvxc,cphi_p,oc,ldim2,ldim2,ndimh)
      cc = transpose(dconjg(cphi))
      call matm(cc,oc,coc,ndimh,ldim2,ndimh)
c--------
c      print *, ' band mode xxx4'
      do iix =1,ndimh  !; print *, ' band mode xxx5 iix=',iix
      do iiy =1,ndimh
c        ovv(iiy,iix) =abs( gpg(iiy,iix)
c     &  +sum( dconjg(cphi(1:ldim2,iiy))*cphi_p(1:ldim2,iix))
c     &   )**2
c        ovv(iiy,iix) =abs( gpg(iiy,iix)
c     &    + sum( dconjg(cphi(1:ldim2,iiy))*
c     &           matmul(ovvx,cphi_p(1:ldim2,iix)) )
c     &                ) !**2
        ovv(iiy,iix) =abs( gpg(iiy,iix) +coc(iiy,iix))
      enddo
            iyiy = maxloc( ovv(1:ndimh,iix),dim=1 )
            write(ifdebug,'(2i3,255i4)') iix ,iyiy,
     &      (int(100*ovv(j,iix)),j=1,ndimh)
      enddo
      deallocate(ppovl2,gpg)

c      print *, ' band mode xxx5'
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxxx1'

      ifill = 0
      do iix = 1,ndimh
!  Denenaracy treatment
        eigdeg = -1d99
        istart = 1
        sumdeg = 0d0
c      print *, ' band mode xxx6'
        do iiy=1,ndimh
          if(iiy==ndimh) then
            ovv(istart:ndimh,iix) = sumdeg + ovv(iiy,iix)
          elseif(abs(eigdeg-evl(iiy))>1d-6.and.iiy/=1) then
            ovv(istart:iiy-1,iix) = sumdeg
            eigdeg = evl(iiy)
            istart = iiy
            sumdeg = ovv(iiy,iix)
          else
            sumdeg = sumdeg + ovv(iiy,iix)
          endif
        enddo
c       print *, ' band mode xxx7'
c            write(ifdebug,*)' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
c            write(ifdebug,'(2i3,255i4)') iix ,iix,
c     &      (int(100*ovv(j,iix)),j=1,ndimh)
! Get the connection ix to iy.
        do
          iyiy = maxloc( ovv(1:ndimh,iix),dim=1 )
          if( ifill(iyiy)==0 ) then
            iiyf(iix) = iyiy
            ifill(iyiy)=1
            ovv(1:iyiy-1    ,iix) = 0d0
            ovv(iyiy+1:ndimh,iix) = 0d0
            exit
          else
            ovv(iyiy,iix) = -1d0
          endif
        enddo
c            write(ifdebug,'(2i3,255i4)') iix ,iyiy,
c     &      (int(100*ovv(j,iix)),j=1,ndimh)
c      print *, ' band mode xxx8'
      enddo
      deallocate(oc,coc,cc)
c      print *, ' end of bndconn',cpusec()
      end

      subroutine matm(a,b,c,n1,n2,n3)
c-- interface for matrix multiplication  c=a*b -------------------------
c assumed size-array
c
c!!! takao 2000 Sep.  If you use this routine with cxml of zgemm(GOTO),
c you might have to allocate a(n1,n2+1) to avoid the stack segmanation fault due to the complilar bug.
c
      integer::  n1,n2,n3
      complex(8) :: a(n1,n2),b(n2,n3), c(n1,n3)
c
c      print *,"sumcheck=",sum(a),sum(b),n1,n2,n3
      call ZGEMM ( "N", "N", n1, n3, n2, (1d0,0d0),
     &             a, n1,
     &             b, n2,
     &             (0d0,0d0), c, n1 )
      end

      subroutine mkppovl2(alat,plat,qlat, ng1,ngvec1,ng2,ngvec2,
     i         nbas, rmax, bas,
     o         ppovl)
C- < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
c
      implicit none
      integer(4) ::  nbas, ng1,ng2,nx(3),
     .        ig1,ig2, ngvec1(3,ng1),ngvec2(3,ng2),
     .         n1x,n2x,n3x,n1m,n2m,n3m
      real(8) :: tripl,rmax(nbas),pi
      real(8) :: plat(3,3),qlat(3,3),
     .  alat,tpibaqlat(3,3),voltot, bas(3,nbas)
C     complex(8) :: img =(0d0,1d0)
      complex(8) :: ppovl(ng1,ng2)
      complex(8),allocatable :: ppox(:,:,:)
c-----------------------------------------------------
C     print *,' mkppovl2:'
      pi        = 4d0*datan(1d0)
      voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      tpibaqlat =  2*pi/alat *qlat
c < G1 | G2 >
      n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
      n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
      n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
      n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
      n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
      n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
c
      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
      ppox = 1d99
      do ig1  = 1, ng1
      do ig2  = 1, ng2
        nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
        if( ppox(nx(1),nx(2),nx(3)) == 1d99 ) then
          call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat,
     i    nx(1:3), ! G2 -G1
     o    ppox( nx(1),nx(2),nx(3)))
        endif
      enddo
      enddo
      do ig1 = 1,ng1
      do ig2 = 1,ng2
        nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
        ppovl(ig1,ig2) = ppox( nx(1),nx(2),nx(3) )
      enddo
      enddo
      deallocate(ppox)
      end
