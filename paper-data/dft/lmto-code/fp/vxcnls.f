      subroutine vxcnls(ri,lcut,nr,np,nlm,nsp,yl,gyl,ylwp,rwgt,wp,
     .  rl,lxcfun,bscal,vl,rep,rmu,buf)
C- Gradient correction to nspher. density on a radial and angular mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   ri    :mesh of points
Ci   lcut  :1 if cutoff exc for small rho to avoid blowup in exc
Ci   nr    :number of radial mesh points
Ci   np    :number of points for angular integration
Ci   nlm   :maximum (l+1)**2
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   yl    :Ylm's tabulated on angular mesh of np points
Ci   gyl   :gradient of yl
Ci   ylwp  :ylwp(ip,ilm) = yl(ip,ilm)*wp(ip)
Ci         :where wp are integration weights for angular mesh
Ci   rwgt  :radial mesh weights
Ci   wp    :angular mesh weights (not needed unless debugging)
Ci   rl    :true density on radial mesh, in Ylm representation
Ci   lxcfun:defines xc functional; see lxcf in structures.h
Ci         :If lxcfun<1000, only nonlocal part of vxc,exc calculated here
Ci   bscal :scale vxc+ - vxc- by bscal (nsp=2 only)
Ci   buf   :Buffer to contain redirected standard output.
Ci         :This routine uses awrite to write output to stdout.
Ci         :It gets redirected to buf if buffer mode is turned on
Co Outputs
Co   vl    :vxcnl is added to vl
Co   rep   :int rho * excnl added to rep
Co   rmu   :int rho * vxcnl added to rmu
Cl Local variables
Cl   lxcfnl: if nonzero, GGA to be evaluated by libxc
Cl         : if zero, GGA is added to an LDA calculated elsewhere
Cl         :          GGA is specified by lxcg
Cl   lxcg  : Specifies GGA, if lxcfnl is zero
Cl         :  0    LSDA
Cl         :  1    Langreth-Mehl
Cl         :  2    PW91
Cl         :  3    PBE
Cl         :  4    PBE with Becke exchange
Cr Remarks
Cr   If vxc+ - vxc- is scaled by bscal, terms rep and rmu cease to
Cr   lose their meaning. the potential is merely scaled; no attempt
Cr   is made to construct a functional corresponding to the scaling.
Cu Updates
Cu   10 Apr 19 Replace write statements with awrite ... buffers parallel
Cu   24 Nov 16 (Jerome Jackson) correction to GGA potential (de/dgradient term included)
Cu             requires additional argument to xc_libxc call
Cu   31 Jul 14 XC B-field can be scaled by bscal
Cu   09 Dec 13 First cut at using libxc functional
Cu   03 Jul 13 Replace f77 pointers with f90 ones
Cu   05 Apr 09 reduced the calling arguments for better portability
Cu   29 Apr 05 (ATP) adaped to new vxcnsp
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lcut,nr,np,nlm,nsp,lxcfun
      double precision bscal
      double precision ri(nr),yl(np,nlm),gyl(np,nlm,3),ylwp(np,nlm),
     .  wp(np),rep(2),rmu(2),rwgt(nr),rl(nr,nlm,nsp),vl(nr,nlm,nsp)
      character*(*) buf
C ... Dynamically allocated local arrays
      real(8), allocatable :: rp(:,:,:),agrp(:,:,:),grp(:,:,:,:),
     .  ggrp(:,:,:),vxcnl(:,:,:),excnl(:,:),agrl(:,:,:)
C ... Local parameters
      double precision pi,vhold,tol,sumnl(0:20),repnl(4),rmunl(4),weight
      double precision rho0(nr,nsp),vnl(nr,nsp),rhor0(2)
      integer xctype(2),lxcg,lxcfnl,lxcf,stdo
      integer ilm,ip,ipr,ir,i,l,ll,lmax,nn,lopgfl,lmx,jx
      integer ir0,ir1,nri,incn
      procedure(integer) :: nglob
      logical lz
      data tol/1d-15/

C     call pshpr(80)
      call tcn('vxcnls')
      call getpr(ipr)
      stdo = nglob('stdo')
      pi = 4d0*datan(1d0)
      nn = 6  ! number of points used to differentiate radial f
      call dpzero(repnl,4)
      call dpzero(rmunl,4)
      allocate(excnl(nr,nlm),vxcnl(nr,nlm,nsp))
      call dpzero(excnl,nr*nlm); call dpzero(vxcnl,nr*nlm*nsp)
      lz = ri(1) == 0
      if (lxcfun >= 1000) then   ! a libxc functional.  Check if need gradients
        call xc_libxc_type(lxcfun,0,xctype)
        if (xctype(1) > 2 .or. xctype(2) > 2)
     .    call rx('vxcnls: functional not implemented')
        if (xctype(1) == 2 .or. xctype(2) == 2) then ! This is GGA
C         lxcfl = 0        ! No potential calculated by vxcnsl (LDA)
          lxcfnl = lxcfun  ! Potential calculated by vxcnls (GGA)
          lxcg = 1         ! Require gradients
          incn = nr
        else ! This is LDA
C         lxcfl = lxcfun ! Evaluate XC functional through vxcnsl (LDA)
          lxcfnl = 0     ! No potential calculated by vxcnsl (GGA)
          lxcg = 0       ! No gradients
          incn = nr
        endif
      else
C       lxcfl = 0               ! Local part is evaluated elsewhere
        lxcfnl = 0              ! No libxc potential calculated by vxcnls
        lxcg = mod(lxcfun/100,100) ! nonlocal part of GGA by vxcnls
        incn = 50
      endif

C --- Make ylwp = yl*wp for fast multiplication ---
C      do  ilm = 1, nlm
C        do  ip = 1, np
C          ylwp(ip,ilm) = yl(ip,ilm)*wp(ip)
C        enddo
C      enddo

C --- Generate density point-wise through sphere ---
      allocate(rp(nr,np,nsp))
      do  i = 1, nsp
        call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0,
     .             rp(1,1,i),nr)
      enddo

C --- If negative density, set to tol ---
      do  i = 1, nsp
        do  ip = 1, np
          do  ir = 1, nr
            if (rp(ir,ip,i) <= 0d0) rp(ir,ip,i) = tol
          enddo
        enddo
      enddo

C --- Potential from spherical part of density (valid for small r) ---
      if (lz) then
        do  i = 1, nsp
          call dpcopy(rl(1,1,i),rho0(1,i),1,nr,dsqrt(4*pi))
          do  ir = 1, nr
            rho0(ir,i) = rho0(ir,i)*ri(ir)**2
          enddo
        enddo
        call dpzero(vnl,nr*nsp)
        if (lxcfun > 1000) then
          lxcf = lxcfun
        else
          lxcf = lxcfun - mod(lxcfun,100)
        endif
        call vxc0sp(lxcf,ri,rwgt,rho0,nr,vnl,excnl,rhor0,1d0,
     .    repnl(3),rmunl(3),nsp)
C        call vxc0gc(nr,nsp,ri,rwgt,rho0,vnl,excnl,repnl(3),rmunl(3),
C     .    100*lxcg)
        call dscal(nr,dsqrt(4*pi),excnl,1)

C        do  i = 1, nsp
C          if (ipr >= 30 .and. i == 1)
C     .      print 4, rmunl(i+2),repnl(i+2),'  (l=0 rho)'
C          if (ipr >= 30 .and. i == 2)
C     .      print 5, rmunl(i+2),repnl(i+2),
C     .      rmunl(3)+rmunl(4),repnl(3)+repnl(4)
C          call dpcopy(vnl(1,i),vxcnl(1,1,i),1,nr,dsqrt(4*pi))
C        enddo
C

        do  i = 1, nsp
          if (ipr >= 30 .and. i == 1)
     .      call awrit2(' vxcnls: nlc rmu=%;11,6D  rep=%;11,6D  (l=0 rho)',buf,100,
     .      stdo,rmunl(i+2),repnl(i+2))

          if (ipr >= 30 .and. i == 2) call awrit4(
     .      ' spin 2:%8f %;11,6D      %;11,6D%N'//
     .      '  total:%8f %;11,6D      %;11,6D',
     .      buf,100,
     .      stdo,rmunl(i+2),repnl(i+2),rmunl(3)+rmunl(4),repnl(3)+repnl(4))
          call dpcopy(vnl(1,i),vxcnl(1,1,i),1,nr,dsqrt(4*pi))
        enddo
C       call prrmsh('l=0 vxcnl in vxcnls',ri,vxcnl,nr,nr,1)
C       call prrmsh('rl in vxcnls',ri,rl,nr,nr,nlm)
C       call prrmsh('wl in vxcnls',ri,rwgt,nr,nr,1)
      endif

C --- Gradient of density point-wise through sphere, and laplacian ---
      lopgfl = 0
      if (lz) lopgfl = 10
      allocate(grp(nr,np,3,nsp),ggrp(nr,np,nsp))
      do  i = 1, nsp
       call gradfl(ll(nlm),nlm,nr,np,1,nr,1,lopgfl,nn,ri,yl,gyl,
     .    rl(1,1,i),grp(1,1,1,i),ggrp(1,1,i))
      enddo

C --- agrl = abs grad rho in its Yl-representation ---
      if (lxcfnl == 0) then
        allocate(agrp(nr,np,nsp),agrl(nr,nlm,nsp))
        do  i = 1, nsp
          do  ip = 1, np
          do  ir = 1, nr
            agrp(ir,ip,i) =
     .      dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
          enddo
          enddo
          call dgemm('N','N',nr,nlm,np,1d0,agrp(1,1,i),nr,ylwp,np,0d0,
     .      agrl(1,1,i),nr)
        enddo
        deallocate(agrp)
      else
        allocate(agrl(1,1,1))
      endif

C --- Do gradient in blocks nr-incn..nr, nr-2*incn..nr-incn ... ---
      ir1  = nr
      lmax = ll(nlm)
   10 continue
        ir0 = max(ir1-incn,1)
        if (lz .and. lxcf < 1000) ir0 = max(ir1-incn,2)
        nri = ir1-ir0+1
        if (nri > 0) then
        vhold = (vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2

C ...   Gradient-corrected Vxc for points between ir0 and ir1
        i = nri*np
        call xxcnls(lxcfnl,lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,ri,
     .    yl,gyl,ylwp,wp,rp,ggrp,grp,agrl,rl,rwgt,lcut,vxcnl,excnl,
     .    sumnl)
        ir1 = ir0-1

C ... Check rmu to determine largest lmax next pass
        do  l = lmax, 0,-1
          lmx = l
          if (dabs(sumnl(l)) > 1d-7) exit
        enddo
        lmax = lmx

        if (dabs(vhold-(vxcnl(ir0,1,1)+vxcnl(ir0,1,nsp))/2)>= 1d-6 .or.
     .    dabs(vhold) <= 1) goto 10
        endif
        deallocate(rp,grp,ggrp,agrl)

C --- Scale nonlocal potential
      if (nsp == 2 .and. bscal /= 1) then
        call bxcscale(nr*nlm,nsp,bscal,nr*nlm,vxcnl)
      endif
C      call prmx('vnl(1)',vxcnl,nr,nr,nlm)
C      call prmx('vnl(2)',vxcnl(1,1,2),nr,nr,nlm)

C --- Nonlocal rho*exc, rho*vxc ---
      do  i = 1, nsp
        do  ilm = 1, nlm
        do  ir = 1, nr
          weight = ri(ir)**2*rwgt(ir)
          rmunl(i) = rmunl(i) + rl(ir,ilm,i)*vxcnl(ir,ilm,i)*weight
          repnl(i) = repnl(i) + rl(ir,ilm,i)*excnl(ir,ilm)*weight
        enddo
        enddo
C        if (ipr >= 30 .and. i == 1) print 4, rmunl(i),repnl(i)
C        if (ipr >= 30 .and. i == 2) print 5, rmunl(i),repnl(i),
C     .    rmunl(1)+rmunl(2),repnl(1)+repnl(2)

        if (ipr >= 30 .and. i == 1)
     .    call awrit2(' vxcnls: nlc rmu=%;11,6D  rep=%;11,6D',buf,100,
     .    stdo,rmunl(i),repnl(i))
        if (ipr >= 30 .and. i == 2) call awrit4(
     .    ' spin 2:%8f %;11,6D      %;11,6D%N'//
     .    '  total:%8f %;11,6D      %;11,6D',
     .    buf,100,
     .    stdo,rmunl(i),repnl(i),rmunl(1)+rmunl(2),repnl(1)+repnl(2))

        rep(i) = rep(i) + repnl(i)
        rmu(i) = rmu(i) + rmunl(i)
      enddo

C --- Add nonlocal vxc into vl ----
      call daxpy(nr*nlm*nsp,1d0,vxcnl,1,vl,1)
      if (lz) then
        do  i = 1, nsp
        vl(1,1,i) = (vl(2,1,i)*ri(3)-vl(3,1,i)*ri(2))/(ri(3)-ri(2))
        jx = 1
        call polint(ri(2),vl(2,1,i),nr-1,nn,ri,0d0,0,jx,vl(1,1,i),vhold)
        if (ipr >= 50 .and. dabs(vhold) > dabs(vl(1,1,i)/100))
     .        print 1,vl(1,1,i),vhold/vl(1,1,i)*100
    1     format(' vxcnls (warning): expect error in V at origin:',
     .  'V=',1pe10.3,' est err=',0pf7.1,'%')
        enddo
      endif
C     call prrmsh('nlocal v(l=0)',ri,vxcnl,nr*nlm,nr,1)

      call tcx('vxcnls')

C --- Print out rmu by angular momentum ---
      if (ipr < 35) return
      lmax = ll(nlm)
      do  i = 1, nsp
        do  l = 0, lmax
          sumnl(l) = 0d0
        enddo
        do  ilm = 1, nlm
          l = ll(ilm)
          do  ir = 1, nr
            sumnl(l) = sumnl(l) + ri(ir)**2*rl(ir,ilm,i)*vxcnl(ir,ilm,i)*rwgt(ir)
          enddo
        enddo
C        if (i == 1) print 2, (sumnl(l),l=0,lmax)
C        if (i == 2) print 3, (sumnl(l),l=0,lmax)
        if (i == 1) call awrit3(' rvnlc by L: %;12,6D%n;10,6D',buf,100,stdo,sumnl,lmax,sumnl(1))
        if (i == 2) call awrit3('     spin 2: %;12,6D%n;10,6D',buf,100,stdo,sumnl,lmax,sumnl(1))
      enddo

C    2 format(' rvnlc by L: ',f12.6,4f10.6:/(15x,4f10.6))
C    3 format('     spin 2: ',f12.6,4f10.6:/(15x,4f10.6))
C    4 format(' vxcnls: nlc rmu=',f11.6,'  rep=',f11.6,a)
C    5 format(' spin 2:         ',f11.6,'      ',f11.6/
C     .       '  total:         ',f11.6,'      ',f11.6)
      deallocate(excnl,vxcnl)
C      if (ipr >= 40) print 887,
C     .  vl(1,1,1), vxcnl(1,1,1), vl(nr,1,1), vxcnl(nr,1,1)
C  887 format(' V_0(0)=',f15.6,'  nloc VXC_0(0)=',f12.6/
C     .       ' V_0(R)=',f15.6,'  nloc VXC_0(R)=',f12.6)

C     call poppr
      end

      subroutine xxcnls(lxcfnl,lxcg,lmax,ir0,ir1,nr,np,nlm,nsp,nn,
     .  ri,yl,gyl,ylwp,wp,rp,ggrp,grp,agrl,rl,rwgt,lcut,vxcnl,excnl,
     .  sumnl)
C- GGA potential, or nonlocal part, inside a sphere for a subset of radii
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcfnl:If > 0 calculate entire XC potential with libxc call
Ci         :If -, calculate nonlocal part of potential according to lxcg
Ci   lxcg  :Used only if lxcfnl=0: calculate nonlocal part of GGA
Ci   lmax  :maximum l for a given site
Ci   ir0   :Evaluate potential for radial mesh points ir0:ir1
Ci   ir1   :Evaluate potential for radial mesh points ir0:ir1
Ci   nr    :number of radial mesh points
Ci   np    :Number of angular mesh points
Ci   nlm   :number of Ylms
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nn    :number of points used to differentiate radial f
Ci         :Used only if lxcfnl=0
Ci   ri    :radial mesh
Ci   yl    :(real) spherical harmonics
Ci   gyl   :gradient of yl
Ci   ylwp  :gl * wp
Ci   wp    :angular mesh weights
Ci   rp    :angular mesh (on 4 pi radians)
Ci   ggrp  :laplacian rho
Ci         :If nsp=2, ggrh(:,1), ggrh(:,2) => laplacians for each spin
Ci   grp   :gradient rho
Ci   agrl  :
Ci   rl    :
Ci   rwgt  :radial mesh weights
Ci   lcut  :
Cl Local variables
Cl   grhop2:|grad rho|^2
Cl         :If spin polarized,|grad rhoup|^2 and |grad rhodn|^2
Cl         :Used only if lxcfnl is nonzero
Cl     ... The following are used only if lxcfnl is zero
Cl   agrp  :grad total rho . grad abs grad total rho
Cl   gagrp :grad abs grad rho -> grad rho . grad abs grad rho
Co Outputs
Co   vxcnl :GGA potential, or nonlocal part if lxcfnl=0
Co   excnl :GGA energy density, or nonlocal part if lxcfnl=0
Co   sumnl
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Nov 16 correction to GGA potential (de/dgradient term included)
Cu             requires additional argument to xc_libxc call
Cu   09 Dec 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ir0,ir1,nr,np,nlm,nsp,nn,lmax,lcut,lxcfnl,lxcg
      double precision ri(nr),rwgt(nr),wp(np),yl(np,nlm),ylwp(np,nlm),
     .  gyl(np,nlm,3),rp(nr,np,nsp),ggrp(nr,np,nsp),
     .  grp(nr,np,3,nsp),vxcnl(nr,nlm,nsp),excnl(nr,nlm),
     .  agrl(nr,nlm,nsp),rl(nr,nlm,nsp),sumnl(0:20)
C ... Dynamically allocated local arrays
      real(8), allocatable :: agrp(:,:,:),gagrp(:,:,:,:),vxcp(:,:,:)
      real(8), allocatable :: excp(:,:),wkl(:,:)
      real(8), allocatable :: exp(:,:),vxp(:,:,:)
      real(8), allocatable :: ecp(:,:),vcp(:,:,:)
      real(8), allocatable :: gradtmp(:,:,:,:),dxdg2(:,:,:)
      real(8), allocatable :: dxdg2yl(:,:),dxtmp(:,:)
C ... Local parameters
      integer :: ll,ir,ip,i,nri,ilm,l,iprint,mlm,lopgfl,dirn
      double precision :: xx(1),tmp
C ... External calls
      external dgemm,dpzero,gradfl,vxcgga,vxnlcc,vxnloc

! --- Setup ---
!     call pshpr(80)
      nri = ir1-ir0+1
      allocate(vxcp(ir0:ir1,np,nsp),excp(ir0:ir1,np))

! --- GGA through libxc call; including ``non-local'' parts
      if (lxcfnl /= 0) then
        allocate(agrp(ir0:ir1,np,2*nsp-1))
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = ir0, ir1
              agrp(ir,ip,i) =
     .        grp(ir,ip,1,i)**2 + grp(ir,ip,2,i)**2 + grp(ir,ip,3,i)**2
            enddo
            if (nsp == 2) then
              do  ir = ir0, ir1
              agrp(ir,ip,3) =
     .          grp(ir,ip,1,1)*grp(ir,ip,1,2) +
     .          grp(ir,ip,2,1)*grp(ir,ip,2,2) +
     .          grp(ir,ip,3,1)*grp(ir,ip,3,2)
              enddo
            endif
          enddo
        enddo

! ...   call libxc for all the points in one blow ... requires ir0 = 1
!       since some arrays dimensioned (1:nr) and others (ir0:ir1)
        if (ir0 /= 1) call rx('xxcnls: fix libxc call for ir0>1')
        allocate(vxp(ir0:ir1,np,nsp),exp(ir0:ir1,np))
        allocate(vcp(ir0:ir1,np,nsp),ecp(ir0:ir1,np))
        allocate(dxdg2(2*nsp-1,nr,np))

        call xc_libxc(nri*np,nsp,-lxcfnl,rp(ir0,1,1),rp(ir0,1,nsp),
     .    xx,xx,xx,agrp(ir0,1,1),xx,agrp(ir0,1,nsp),
     .    agrp(ir0,1,2*nsp-1),xx,xx,xx,exp,ecp,vxp,vcp,
     .    vxp(ir0,1,nsp),vcp(ir0,1,nsp),
     .    excp(ir0,1),vxcp(ir0,1,1),vxcp(ir0,1,nsp),dxdg2(1,1,1))
CJJ, require d(exc)/d(gradient rho^2) from libxc call
C see, eg, equation (14) of libxc paper
C Marques, Oliveira, Burnus, CPC 183 (2012) 2272--2281
C
C dxdg2 corresponds to the libxc variable vsigma
C
        allocate(dxtmp(nr,np),gradtmp(nr,np,3,3),dxdg2yl(nr,nlm)) ! require dxdg also in ylm repn for gradfl function
        do i = 1, 2*nsp-1
          dxtmp(:,:)=dxdg2(i,:,:) !awkward array indexing from libxc
          call dgemm('N','N',nr,nlm,np,1.d0,dxtmp,nr,ylwp,np,0.d0,dxdg2yl,nr) ! dxdg2 in (r,p), need ylwp here ;-)
          ! grad(de/dsigma_i)
          call gradfl(lmax,nlm,nr,np,ir0,ir1,0,10,nn,ri,yl,gyl,dxdg2yl,gradtmp(1,1,1,i),0.d0) ! gradtmp in (r,p)
        end do

        !2*grad(de/dsigma_uu).grad(rho_u) + grad(de/dsigma_ud).grad(rho_d) + 2*(de/dsigma_uu)*lap(rho_u) + (de/dsigma_ud)*lap(rho_d), etc
        do ir=1,nr
          do ip=1,np
            ! easiest to consider spins separately
            ! spin 1
            tmp=0.d0
            do dirn=1,3
               tmp=tmp+2*gradtmp(ir,ip,dirn,1)*grp(ir,ip,dirn,1) ! grad and lap rho calculated previously
               if (nsp==2) tmp=tmp+gradtmp(ir,ip,dirn,2)*grp(ir,ip,dirn,2)
            end do
            tmp=tmp+2*dxdg2(1,ir,ip)*ggrp(ir,ip,1)
            if(nsp==2)tmp=tmp+dxdg2(2,ir,ip)*ggrp(ir,ip,2)
            vxcp(ir,ip,1)=vxcp(ir,ip,1)-tmp

            ! spin 2
            if(nsp==2)then
              tmp=0.d0
              do dirn=1,3
                 tmp=tmp+2*gradtmp(ir,ip,dirn,3)*grp(ir,ip,dirn,2)+gradtmp(ir,ip,dirn,2)*grp(ir,ip,dirn,1)
              end do
              tmp=tmp+2*dxdg2(3,ir,ip)*ggrp(ir,ip,2)+dxdg2(2,ir,ip)*ggrp(ir,ip,1)
              vxcp(ir,ip,2)=vxcp(ir,ip,2)-tmp
            end if ! nsp==2
          end do
        end do
        deallocate(gradtmp,dxdg2,dxdg2yl,dxtmp)
CJJ
        deallocate(agrp,exp,vxp,ecp,vcp)

! --- nonlocal part of GGA from vxcgga or vxcnloc
      else
        allocate(agrp(ir0:ir1,np,3*nsp-2),gagrp(ir0:ir1,np,3,nsp))

! --- gagrp(store in vxcp) = grad rho . grad abs grad rho ---
        lopgfl = 0
        if (ri(ir0) == 0) lopgfl = 10
        do  i = 1, nsp
          call gradfl(lmax,nlm,nr,np,ir0,ir1,0,lopgfl,nn,ri,yl,gyl,
     .      agrl(1,1,i),gagrp(ir0,1,1,i),0d0)
        enddo
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = ir0, ir1
              vxcp(ir,ip,i) = gagrp(ir,ip,1,i)*grp(ir,ip,1,i)
     .                       +gagrp(ir,ip,2,i)*grp(ir,ip,2,i)
     .                       +gagrp(ir,ip,3,i)*grp(ir,ip,3,i)
            enddo
          enddo
        enddo

!   --- store in agrp:  grad total rho . grad abs grad total rho ---
        if (nsp == 2) then
          do  ip = 1, np
          do  ir = ir0, ir1
            agrp(ir,ip,1) =
     .        (grp(ir,ip,1,1)+grp(ir,ip,1,2))*
     .        (gagrp(ir,ip,1,1)+gagrp(ir,ip,1,2)) +
     .        (grp(ir,ip,2,1)+grp(ir,ip,2,2))*
     .        (gagrp(ir,ip,2,1)+gagrp(ir,ip,2,2)) +
     .        (grp(ir,ip,3,1)+grp(ir,ip,3,2))*
     .        (gagrp(ir,ip,3,1)+gagrp(ir,ip,3,2))
          enddo
          enddo
        endif

!   --- Copy grad rho . grad abs grad rho into gagrp ---
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = ir0, ir1
              gagrp(ir,ip,i,1) = vxcp(ir,ip,i)
            enddo
          enddo
        enddo
        if (nsp == 2) then
          do  ip = 1, np
            do  ir = ir0, ir1
              gagrp(ir,ip,3,1) = agrp(ir,ip,1)
            enddo
          enddo
        endif
C       call px('gr.gagr',nri,nlm,1,np,ri(ir0),wp,gagrp(ir0,1,3,1),yl,wkl)

!   --- Make agrp+,agrp- for ir0 .. ir1 ---
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = ir0, ir1
            agrp(ir,ip,i) =
     .      dsqrt(grp(ir,ip,1,i)**2+grp(ir,ip,2,i)**2+grp(ir,ip,3,i)**2)
            enddo
          enddo
        enddo

!   --- Make agrp (total rho) agrp+.agrp-  for ir0 .. ir1 ---
        if (nsp == 2) then
          do  ip = 1, np
            do  ir = ir0, ir1
              agrp(ir,ip,3) =
     .          dsqrt((grp(ir,ip,1,1)+grp(ir,ip,1,2))**2 +
     .                (grp(ir,ip,2,1)+grp(ir,ip,2,2))**2 +
     .                (grp(ir,ip,3,1)+grp(ir,ip,3,2))**2)
              agrp(ir,ip,4) =
     .                 grp(ir,ip,1,1)*grp(ir,ip,1,2) +
     .                 grp(ir,ip,2,1)*grp(ir,ip,2,2) +
     .                 grp(ir,ip,3,1)*grp(ir,ip,3,2)
            enddo
          enddo
C       call px('x',nri,nlm,nsp,np,ri(ir0),wp,agrp(ir0,1,3),yl,wkl)
        endif

C --- Make nonlocal potential for points ir0 .. ir1 ---
        call dpzero(excp,nri*np); call dpzero(vxcp,nri*np*nsp)
        do  ip = 1, np
          i = 1 ; if (nsp == 2) i = 4 ! Enables vxcnsp to pass bounds check
          if (lxcg > 2) then
            call vxcgga(lxcg,nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),
     .        agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1),
     .        ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,i),
     .        gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1),
     .        gagrp(ir0,ip,nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp),
     .        excp(ir0,ip))
          elseif (lcut == 0) then
            call vxnloc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),
     .        agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1),
     .        ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,i),
     .        gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1),gagrp(ir0,ip,
     .        nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp),excp(ir0,ip))
          else
            call vxnlcc(nri,nsp,rp(ir0,ip,1),rp(ir0,ip,nsp),
     .        agrp(ir0,ip,1),agrp(ir0,ip,nsp),ggrp(ir0,ip,1),
     .        ggrp(ir0,ip,nsp),agrp(ir0,ip,2*nsp-1),agrp(ir0,ip,i),
     .        gagrp(ir0,ip,2*nsp-1,1),gagrp(ir0,ip,1,1),gagrp(ir0,ip,
     .        nsp,1),vxcp(ir0,ip,1),vxcp(ir0,ip,nsp),excp(ir0,ip))
          endif
        enddo
        deallocate(agrp,gagrp)
      endif  ! generate of potential pointwise through sphere

C     call bxcscale(nri*np,nsp,bscal,nri*np,vxcp)

C ... (test): yl projection of various quantities'
C      call px('rho',nr,nlm,nsp,np,ri,wp,rp,yl,wkl)
C      call px('ggrh',nr,nlm,nsp,np,ri,wp,ggrp,yl,wkl)
C      call px('agrh',nri,nlm,nsp,np,ri(ir0),wp,agrp,yl,wkl)
C      call px('gr.gagr',nri,nlm,nsp,np,ri(ir0),wp,gagrp,yl,wkl)
C      call px('vxc',nri,nlm,nsp,np,ri(ir0),wp,vxcp,yl,wkl)
C      call px('exc',nri,nlm,nsp,np,ri(ir0),wp,excp,yl,wkl)

C --- Yl-projection of vxcp,excp into vxcnl,excnl ---
      allocate(wkl(ir0:ir1,nlm))
      mlm = (lmax+1)**2
      call dgemm('N','N',nri,mlm,np,1d0,excp(ir0,1),nri,ylwp,np,0d0,
     .  wkl,nri)
      do  ilm = 1, mlm
        do  ir = ir0, ir1
          excnl(ir,ilm) = wkl(ir,ilm)
        enddo
      enddo
      do  i = 1, nsp
        call dgemm('N','N',nri,mlm,np,1d0,vxcp(ir0,1,i),nri,ylwp,np,0d0,
     .    wkl,nri)
        do  ilm = 1, mlm
          do  ir = ir0, ir1
            vxcnl(ir,ilm,i) = wkl(ir,ilm)
          enddo
        enddo
      enddo
C     call prmr(nr,ri,vxcnl,nlm)

C --- Estimate rmu in this ir0 ir1 interval by angular momentum ---
      do  i = 1, nsp
        do  l = 0, lmax
          sumnl(l) = 0d0
        enddo
        do  ilm = 1, nlm
          l = ll(ilm)
          do  ir = ir0, ir1
            sumnl(l) = sumnl(l) +
     .        rl(ir,ilm,i)*vxcnl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
          enddo
        enddo
        if (iprint() >= 80) then
          if (i == 1) print 1,ri(ir0),(sumnl(l),l = 0,lmax)
          if (i == 2) print 2,(sumnl(l),l = 0,lmax)
    1     format(' R>',f8.6,': ',f12.6,4F10.6:/(15x,4F10.6))
    2     format('     spin 2: ',f12.6,4F10.6:/(15x,4F10.6))
        endif
      enddo
C     call poppr

      deallocate(vxcp,excp,wkl)

      end
C      subroutine px(strn,nr,nlm,nsp,np,ri,wp,fp,yl,fl)
C      implicit none
C      character *(*) strn
C      integer nr,nlm,nsp,np
C      double precision ri(*),wp(np),fp(nr,np,nsp),yl(np,nlm),
C     .  fl(nr,nlm,nsp)
C
C      call dpzero(fl,nr*nlm*nsp)
C      call fp2yl(nr,nlm,nsp,np,wp,fp,yl,0d0,fl)
C      print *, fl(1,1,1), fl(1,1,2)
C      print *, fl(591,1,1), fl(591,1,2)
C      print *, strn
C      call prmr(nr,ri,fl,nlm)
C
C      end
C
C      subroutine prmr(nr,ri,f,nl)
C      implicit none
C      integer nr,nl,ir,j,fopna,ifi,ir0
C      double precision ri(nr),f(nr,nl)
C      character*(10) fmt
C      ifi = fopna('out',19,0)
C      ir0 = 1
CC ... first nonzero l=0 point ...
C      do  20  ir = 1, nr
C        ir0 = ir
C        if (dabs(f(ir,1)) > 1d-12) goto 21
C   20 continue
C   21 continue
C
C      write(ifi,*) nr-ir0+1, nl+1
C      do  10  ir = ir0, nr
CC        write(ifi,333) ri(ir), (f(ir,3*j-2), j=1, nl)
C        write(ifi,333) ri(ir), (f(ir,j), j=1, nl)
CC 333   format(f12.7,(7g18.10:/12x))
C  333   format(f12.9,(9g18.10:/12x))
C   10 continue
C      call fclose(ifi)
C      print *, 'prmr:'
C      pause
C      end
C
