      subroutine vxc0sp(lxcfun,rofi,rwt,rho,nr,v,exc,rho0,bscal,rep,rmu,
     .  nsp)
C- Adds xc part to spherical potential, makes integrals rmu and rep
C ----------------------------------------------------------------------
Ci Inputs
Ci   rofi  :radial mesh points
Ci   rwt   :radial mesh weights
Ci   rho   :density = (true rho)*(4*pi*r**2)
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   bscal :factor to scale effective xc magnetic field (for const.B)
Co Outputs
Co   v     :vxc is added to v
Co   exc   :XC energy density is returned
Co   rho0  :density extrapolated to origin
Co   rep   :integral rho * exc.
Co   rmu   :integral rho * vxc.
Cl Local variables
Cl  lxcf     :type xc functional; see structures.h
Cr Remarks
Cu Updates
Cu   23 Nov 16 (Jerome Jackson) added d(e_xc)/d(grad rho) term to libxc GGA potential
Cu   08 Dec 13 Enable calls to libxc library. Argument list was modified.
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   25 Apr 12 (K Belashenko) new bscal; mesh parameter b eliminated
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   18 Jun 04 lxcfun is no longer used and should be deleted from cmd-line
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,lxcfun
      double precision rofi(nr),rwt(nr),v(nr,nsp),rho(nr,nsp),bscal,
     .  rep(nsp),rmu(nsp),rho0(2),qs(2),exc(nr)
C ... Dynamically allocated local arrays
      real(8), allocatable :: ex(:),ec(:),vx(:),vc(:),vx2(:),vc2(:)
      real(8), allocatable :: grh(:,:),ggrh(:),agrh(:,:),grgag(:)
C ... Local parameters
      integer iprint,nglob,stdo
      integer i,ir,isp,lxcf,lxcg,xctype(2)
      double precision pi,rho2,rho3,ub4pi,xx(1),
     .  vxc(nr,2),rp(nr,2),repl(2),rmul(2)
      character*2 st
      double precision, allocatable ::dxdg2(:,:), tmp(:), gradtmp(:)

      stdo = nglob('stdo')
      lxcf = lxcfun
      if (lxcf >= 1000) then   ! a libxc functional.  Check if need gradients
        call xc_libxc_type(lxcf,0,xctype)
        if (xctype(1) > 2 .or. xctype(2) > 2)
     .    call rx('vxc0sp: functional not implemented')
        lxcg = 0
        if (xctype(1) == 2 .or. xctype(2) == 2) lxcg = 1 ! Require gradients
      else
        lxcg = lxcf/100
        lxcf = mod(lxcf,100)
      endif
      pi = 4d0*datan(1d0)
      ub4pi = 1d0/(4d0*pi)
C     intopt = 10*nglob('lrquad')

C --- Add background rho to calculate vxc ---
*      rhobg = 0d0
*      call getsyv('rhobg',rhobg,i)
*      call addzbk(rofi,nr,1,nsp,rho,rhobg,1d0)

C --- Extrapolate rho to origin ---
      do  isp = 1, nsp
      rep(isp) = 0d0
      rmu(isp) = 0d0
      rho2 = rho(2,isp)/rofi(2)**2
      rho3 = rho(3,isp)/rofi(3)**2
      rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
      enddo

C --- Make true rho ---
      do  isp = 1, nsp
      rp(1,isp) = rho0(isp)
        do  ir = 2, nr
          rp(ir,isp) = rho(ir,isp)*ub4pi/rofi(ir)**2
        enddo
C       call ratint(rofi(2),rp(2,isp),4,rofi(1),rp(1,isp),xx)
      enddo

C --- Generate vxc,exc on a mesh ---
      allocate(ex(nr),ec(nr),vx(nr),vc(nr),vx2(nr),vc2(nr))

      if (lxcf > 1000) then ! A libxc functional
        if (lxcg == 0) then ! LDA; no gradients needed
          call xc_libxc(nr,nsp,-lxcf,rp,rp(1,nsp),xx,xx,xx,xx,xx,xx,xx,
     .      xx,xx,xx,ex,ec,vx,vc,vx2,vc2,exc,vxc,vxc(1,nsp),xx)
C         print *, sngl(exc(100)),sngl(vxc(100,1)),sngl(vxc(100,2))
        else
          allocate(dxdg2(2*nsp-1,nr)) ! funky ordering for libxc's interface
          allocate(grh(nr,nsp),agrh(nr,2*nsp-1))
          !vxc=0.d0 ! -JJ, this doesn't seem to be initialised anywhere
          call vxcgr2(-1,nr,nsp,nr,rofi,rp,grh,xx,agrh,xx,xx,xx)
          call xc_libxc(nr,nsp,-lxcf,rp,rp(1,nsp),xx,xx,xx,agrh,
     .      xx,agrh(1,nsp),agrh(1,2*nsp-1),xx,xx,xx,ex,ec,vx,vc,vx2,
     .      vc2,exc,vxc,vxc(1,nsp),dxdg2)

CJJ, inclusion of \vec{grad}.( d(e_xc)/d(grad rho)) ) contribution to V
C see, eg, equation (14) of libxc paper
C Marques, Oliveira, Burnus, CPC 183 (2012) 2272--2281
C
C dxdg2 == libxc's vsigma == de/(grad(rho)^2), use chain rule
C
C only radial dependence: evaluated in spherical coords
          allocate(tmp(nr),gradtmp(nr))
          do i=1,nsp
            ! note nspin=2 case mixes spins
            ! tmp =(de/dsigma) * dsigma/d(grad rho); use chain rule
            if(i==1)then
              tmp(:)=2*dxdg2(1,:)*grh(:,1)
              if(nsp==2) then
                tmp(:)=tmp(:) + dxdg2(2,:)*grh(:,2)
              end if
            else
              tmp(:)=2*dxdg2(3,:)*grh(:,2) + dxdg2(2,:)*grh(:,1)
            end if
            call radgrx(nr,nr,1,rofi,tmp,gradtmp)
            do ir=2,nr ! note 1/r: first value diverges
                vxc(ir,i)=vxc(ir,i)-(gradtmp(ir)+2*tmp(ir)/rofi(ir))
            end do
            ! extrapolate...
            vxc(1,i)=vxc(1,i)- gradtmp(2)+2*tmp(2)/rofi(2)
     .        -(gradtmp(3)+2*tmp(3)/rofi(3)-gradtmp(2)-2*tmp(2)/rofi(2))
     .        *rofi(2)/(rofi(3)-rofi(2))
            ! or, leave vxc(origin) uncorrected by this term
          end do
          deallocate(tmp,gradtmp,dxdg2,grh,agrh)
        endif

        lxcg = 0                ! No further gradient corrections

      elseif (lxcf > 2) then
        call evxcp(rp,rp(1,2),nr,nsp,lxcf,ex,ec,exc,vx,vx2,vc,vc2,vxc,
     .    vxc(1,2))
        do  isp = 1, nsp
          vxc(1,isp) = (vxc(2,isp)*rofi(3)-vxc(3,isp)*rofi(2))
     .                 /(rofi(3)-rofi(2))
        enddo
      elseif (lxcf > 0) then
        if (nsp == 1) then
          call evxcv(rp,rp,nr,nsp,lxcf,exc,ex,ec,vxc,vx,vc)
        else
C          call xc_libxc(nr,nsp,-65553,rp,rp(1,nsp),xx,xx,xx,xx,xx,xx,xx,
C     .      xx,xx,xx,ex,ec,vx,vc,vx2,vc2,exc,vxc,vxc(1,nsp))
          call dpadd(rp(1,2),rp,1,nr,1d0)
          call evxcv(rp(1,2),rp,nr,2,lxcf,exc,ex,ec,vxc,vx,vc)
          call dpadd(rp(1,2),rp,1,nr,-1d0)
          call dpadd(rp,rp(1,2),1,nr,1d0)
          call evxcv(rp,rp(1,2),nr,2,lxcf,exc,ex,ec,vxc(1,2),vx,vc)
          call dpadd(rp,rp(1,2),1,nr,-1d0)
        endif
      else
        call dpzero(exc,nr)
        call dpzero(vxc,nr*2)
      endif
      deallocate(ex,ec,vx,vc,vx2,vc2)

C --- Add gradient correction ---
      if (lxcg /= 0) then
        do  i = 1, nsp
          repl(i) = 0d0; rmul(i) = 0d0
          do  ir = 1, nr
            repl(i) = repl(i) + rwt(ir)*rho(ir,i)*exc(ir)
            rmul(i) = rmul(i) + rwt(ir)*rho(ir,i)*vxc(ir,i)
          enddo
        enddo
        allocate(grh(nr,nsp),ggrh(nr*nsp))
        allocate(agrh(nr,3*nsp-2),grgag(nr*(2*nsp-1)))
        call vxcgr2(lxcg,nr,nsp,nr,rofi,rp,grh,ggrh,agrh,grgag,exc,
     .    vxc)
        deallocate(grh,ggrh,agrh,grgag)
      endif

C ... Scale XC field
      if (bscal /= 1d0 .and. nsp == 2) then
C       print *, vxc(nr,1)+vxc(nr,2), vxc(nr,1)-vxc(nr,2)
        call bxcscale(nr,nsp,bscal,nr,vxc)
C       print *, vxc(nr,1)+vxc(nr,2), vxc(nr,1)-vxc(nr,2)
      endif

C --- Integrals ---
      do  i = 1, nsp
        qs(i)  = 0d0
        rep(i) = 0d0
        rmu(i) = 0d0
        do  ir = 1, nr
          qs(i)  = qs(i)  + rwt(ir)*rho(ir,i)
          rep(i) = rep(i) + rwt(ir)*rho(ir,i)*exc(ir)
          rmu(i) = rmu(i) + rwt(ir)*rho(ir,i)*vxc(ir,i)
        enddo
      enddo

C --- Add to V ---
      call dpadd(v,vxc,1,nr,1d0)
      if (nsp == 2) call dpadd(v(1,2),vxc(1,2),1,nr,1d0)

C --- Undo background rho for purposes of calculating vxc ---
*     call addzbk(rofi,nr,1,nsp,rho,rhobg,-1d0)

      if (iprint() < 80) return
      if (lxcg == 0) write(stdo,1)
      if (lxcg /= 0) write(stdo,2)
    1 format(/' vxc0sp: reps(l)     rmu(l)')
    2 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
      do  i = 1, nsp
      st = ' '
      if (i < nsp) st = 'up'
      if (i == 2)   st = 'dn'
      if (lxcg == 0) write(stdo,3) st, rep(i),  rmu(i)
      if (lxcg /= 0) write(stdo,3) st, repl(i), rmul(i),
     .  rep(i)-repl(i), rmu(i)-rmul(i)
      enddo
      if (nsp == 2 .and. lxcg == 0)
     .  write(stdo,3) '  ', rep(1)+rep(2), rmu(1)+rmu(2)
      if (nsp == 2 .and. lxcg /= 0)
     .  write(stdo,3) '  ', repl(1)+repl(2), rmul(1)+rmul(2),
     .  rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)
    3 format(1x,a2,2x,4f12.6)

      end

      subroutine bxcscale(np,nsp,bscal,ndv,vxc)
C- Scale vxc+ - vxc- by bscal
C ----------------------------------------------------------------------
Ci Inputs
Ci   np    : Number of points
Ci   nsp   : vxcscale only operates in nsp=2
Ci   bscal : scale factor
Ci   ndv   : Dimensions vxc
Cio Inputs/Outputs
Cio  vxc   : potential.
Cio        : vxc+ + vxc- is not affected
Cio        : vxc+ - vxc- scaled by bscal
Cr Remarks
Cr
Cu Updates
Cu   24 Jul 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer np,nsp,ndv
      double precision bscal,vxc(ndv,nsp)
C ... Local parameters
      real(8), allocatable :: vwk(:)

      if (bscal == 1 .or. nsp /= 2) return

      allocate(vwk(np))
      call dpcopy(vxc,vwk,1,np,.5d0)
      call dpadd(vwk,vxc(1,2),1,np,-0.5d0) ! bxc
      call dpadd(vxc,vwk,1,np,bscal-1)
      call dpadd(vxc(1,2),vwk,1,np,1-bscal)
      deallocate(vwk)

      end

C      subroutine bxcscalex(np,nsp,bscal,ndv,vxc)
CC- Scale vxc+ - vxc- by bscal
CC ----------------------------------------------------------------------
CCi Inputs
CCi   np    : Number of points
CCi   nsp   : vxcscale only operates in nsp=2
CCi   bscal : scale factor
CCi   ndv   : Dimensions vxc
CCio Inputs/Outputs
CCio  vxc   : potential.
CCio        : vxc+ + vxc- is not affected
CCio        : vxc+ - vxc- scaled by bscal
CCr Remarks
CCr
CCu Updates
CCu   24 Jul 14 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer np,nsp,ndv
C      double precision bscal,vxc(ndv,nsp)
CC ... Local parameters
C      real(8), allocatable :: vwk(:)
C
C      if (bscal == 1 .or. nsp /= 2) return
C
C      allocate(vwk(np))
C      call dpcopy(vxc,vwk,1,np,.5d0)
C      call dpadd(vwk,vxc(1,2),1,np,-0.5d0) ! bxc
C      call dpadd(vxc,vwk,1,np,bscal-1)
C      call dpadd(vxc(1,2),vwk,1,np,1-bscal)
C      deallocate(vwk)
C
C      end
