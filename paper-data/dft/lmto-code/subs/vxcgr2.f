      subroutine vxcgr2(lxcg,nr,nsp,nrx,rofi,rp,
     .  grh,ggrh,agrh,grgagr,exc,vxc)
C- Radial gradients of spherical density and optionally some GGA
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcg  :>0 Make gradients for, call vxcgga
Ci         :-1 Make gradients for call libxc GGA
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nrx   :leading dimension of the radial function arrays
Ci   rofi  :radial mesh points
Ci   rp    :density rho on a radial mesh
Co Outputs
Co   grh   :radial grad rho
Co         :If nsp=2, grh(:,1), grh(:,2) => gradients for each spin
Co  ... Other gradients depend on lxcg.
Co         :Case lxcg=-1, nsp=1
Co   agrh  :|grad rho|^2 = |grh|^2 if nsp=1
C          :If nsp=2
Co         :agrh(:,1) = |grh(:,1)|^2  agrh(:,2)= |grh(:,2)|^2
Co         :agrh(:,3)=  grh(:,1).grh(:,2)
Co         :Case lxcg>0
Co   ggrh  :laplacian rho, radial part
Co         :If nsp=2, ggrh(:,1), ggrh(:,2) => laplacians for each spin
Co   agrh  :abs(grh)
Co         :If nsp=2,
Co         : agrh(:,1)= |grh(:,1)|  agrh(:,2)= |grh(:,2)|
Co         : agrh(:,3)= |grh(:,1)+grh(:,2)|
Co         : agrh(:,4)= grh(:,1)*grh(:,2)
Ci   grgagr:grad rho . grad abs grad rho
Co   exc   :gradient contribution to energy added to exc
Co   vxc   :gradient contribution to potential added to vxc
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   01 Dec 13 Add lxcg to argument list
Cu   18 Jun 04 Bug fix
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,nrx,lxcg
      double precision rp(nrx,nsp),grh(nrx,2),ggrh(nrx,2),
     .  agrh(nrx,4),grgagr(nrx,3),exc(nrx),vxc(nrx,2),rofi(nr)
C ... Local parameters
      integer ir,i
C ... External calls
      external poldvm,polint,radgrx,rx,vxcgga,vxnloc

C     call prrmsh('rho',rofi,rp,nrx,nr,nsp)

C --- grad(rho) ---
      call radgrx(nr,nrx,nsp,rofi,rp,grh)
C     call prrmsh('grh',rofi,grh,nrx,nr,nsp)

C --- Laplacian rho ---
      if (lxcg > 0) then
        call radgrx(nr,nrx,nsp,rofi,grh,ggrh)
C       First two points doubtful ... should extrapolate
C       call prrmsh('ggrh',rofi,ggrh,nrx,nr,nsp)
C       The term below dominates anyway, so no matter
        do  i = 1, nsp
          do  ir = 2, nr
            ggrh(ir,i) = ggrh(ir,i) + 2d0*grh(ir,i)/rofi(ir)
          enddo
          ggrh(1,i) =
     .      (rofi(3)*ggrh(2,i)-rofi(2)*ggrh(3,i))/(rofi(3)-rofi(2))
C         call prrmsh('ggrh',rofi,ggrh,nrx,nr,nsp)

C     --- |grad rho|, |grad rho| . grad |grad rho| ---
          if (lxcg > 0) then
            do  ir = 1, nr
              agrh(ir,i) = dabs(grh(ir,i))
            enddo
            call radgrx(nr,nrx,1,rofi,agrh(1,i),grgagr(1,i))
            do  ir = 1, nr
              grgagr(ir,i) = grh(ir,i)*grgagr(ir,i)
            enddo
          endif
        enddo
      endif

C --- Extra terms g(n), g(n+).g(n-), g(n).g(abs(g(n))) if spin pol ---
      if (nsp == 2 .and. lxcg > 0) then
        do  ir = 1, nr
          agrh(ir,3) = dabs(grh(ir,1)+grh(ir,2))
        enddo
        call radgrx(nr,nrx,1,rofi,agrh(1,3),grgagr(1,3))
        do  ir = 1, nr
          grgagr(ir,3) = (grh(ir,1)+grh(ir,2))*grgagr(ir,3)
        enddo
        do  ir = 1, nr
          agrh(ir,4) = grh(ir,1)*grh(ir,2)
        enddo
      elseif (nsp == 2) then  ! Gradients for libxc, nsp=2
        do  ir = 1, nr
          agrh(ir,1) = (grh(ir,1))**2
          agrh(ir,2) = (grh(ir,2))**2
          agrh(ir,3) = grh(ir,1)*grh(ir,2)
        enddo
        return
      elseif (lxcg < 0) then ! Gradients for libxc, nsp=1
        do  ir = 1, nr
          agrh(ir,1) = (grh(ir,1))**2
        enddo
        return
      endif

C --- Gradient term for all points ---
      if (lxcg >= 3) then
C        lxcf = mod(nglob('lxcf'),100)
C        if (lxcf /= 3 .and. lxcf /= 4) call
C     .    rx('vxcgf2: inconsistent use of local and GGA functionals')
        call vxcgga(lxcg,nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp),
     .    ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4),
     .    grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp),
     .    vxc(1,1),vxc(1,nsp),exc)
      elseif (lxcg == 2) then
        call rx('PW91 no longer implemented')
      else
        call vxnloc(nr,nsp,rp,rp(1,nsp),agrh,agrh(1,nsp),
     .    ggrh,ggrh(1,nsp),agrh(1,2*nsp-1),agrh(1,4),
     .    grgagr(1,2*nsp-1),grgagr,grgagr(1,nsp),
     .    vxc(1,1),vxc(1,nsp),exc)
      endif
      do  i = 1, nsp
        vxc(1,i) = (vxc(2,i)*rofi(3)-vxc(3,i)*rofi(2))/(rofi(3)-rofi(2))
      enddo

      end
      subroutine radgrx(nr,nrx,nsp,ri,f,gf)
      implicit none
      integer nr,nrx,nsp,nn,i,iprint,jx
      double precision ri(nr),f(nrx,nsp),gf(nrx,nsp),tol,egf0
      logical lerr
      parameter (tol=1d-12,nn=6)

      do  i = 1, nsp
C       call prmr(nr,ri,gf,1)
        call poldvm(ri(2),f(2,i),nr-1,nn,.false.,tol,lerr,gf(2,i))
        jx = 1
        call polint(ri(2),gf(2,i),nr-1,nn,ri,0d0,0,jx,gf(1,i),egf0)
        if (iprint() >= 50 .and. dabs(egf0/gf(1,i)) > 1d-2)
     .    print 1, gf(1,i), egf0/gf(1,i)*100
    1   format(' radgrx (warning): expect error in gradient at origin:',
     .    'f=',1pe10.3,' est err=',0pf7.1,'%')
      enddo
      end
