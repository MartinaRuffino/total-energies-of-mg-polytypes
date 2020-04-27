      subroutine solhsg(ph,rsm,eh,lmxh,iopt,cg0,hkl,ghkl)
C- Solid sm Hankels H_kL, or related function, gradient and lapl, at point
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :coordinates relative to centre of function hkl
Ci   rsm   :vector of l-dependent smoothing radii of the function
Ci         :rsm<0 => rsm, eh are independent of l.  Use |rsm(0)| for all l
Ci         :rsm>0 => rsm,eh must be specified for 0..lmax
Ci   eh    :l-dependent Hankel energies
Ci   lmxh  :max angular momentum
Ci   iopt  :Determines the function to return; see Remarks
Ci         :0 return hkl
Ci         :1 return energy derivative of hkl
Ci         :2 return hkl + cg0 * gkl
Ci         :3 return hkl - hkl(rsm=0)
Ci         :4 return hkl - hkl(rsm=0) + cg0 * gkl
Ci   cg0   :coefficients to amout of gkl to add to hkl; see Remarks
Co Outputs
Co   hkl   : H_kL(X), k = 0,1, L = 1,...,(lmax+1)^2
Co   ghkl  :\grad H_L(X), L = 1,...,(lmax+1)^2
Cl Local variables
Cl  fixdrs :T => rsm are l-dependent
Cr Remarks
Cr   Functions are designed for optimally chosen screened basis set.
Cr     Use iopt 0,1  for sm Hankel Hkl or energy derivative
Cr              2-4  for hkL + cg0*gkL or similar.
Cr                   Note that cg0, rsm can be by an automatic
Cr                   prescription; see rsfit.f
Cr
Cr   Gradients of H_kl require H_k+1,l, and H_k,l-1,H_k,l+1
Cb Bugs
Cr   Routine could be greatly simplified by patterning it after solhjg
Cu Updates
Cu   26 Aug 11 Additional kinds of functions may be calculated
Cu   20 Apr 11 Adapted from solgsg
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer lmxh,iopt
      double precision ph(3),rsm(0:lmxh),eh(0:lmxh),cg0(0:lmxh)
C Output parameters
      double precision hkl(0:1,(lmxh+1)**2),ghkl(3,(lmxh+1)**2)
C Local variables
      logical fixdrs,ldot
      integer lmax,isw,l,ilm,ll,ilr,jhd,nlmp
      double precision yl((lmxh+2)**2),ddot,px(3),xx
      double precision tol
      parameter (tol=1d-14)
      double precision rsm0,rsx,r1,ehx,eh0
      double precision hs(0:lmxh+2),dhs(0:lmxh+2),ddhs(0:lmxh+2)
      double precision hsp(0:lmxh+2),dhsp(0:lmxh+2),ddhsp(0:lmxh+2)
      double precision hloc(0:1,(lmxh+2)**2)
      double precision gloc(0:1,(lmxh+2)**2),ggrad(3,(lmxh+1)**2)
      integer kz,kx1,kx2,ky1,ky2
      double precision cz,cx1,cx2,cy1,cy2
C     double precision gkl(0:1,0:lmxh+2)

      ldot = iopt == 1

C ... Spherical harmonics (unit radius)
      r1 = dsqrt(ddot(3,ph,1,ph,1))
      if (r1 > 0) call dpscop(ph,px,3,1,1,1/r1)
      call ropyln(1,px(1),px(2),px(3),lmxh+1,1,yl,xx)

C ... if r1 = 0, only l=0 term survives in H_kL and
C                     l=1 terms survive in \grad G_L
      lmax = lmxh
      if (dabs(r1) < tol .and. lmxh > 1) then
        lmax = 1
        call dpzero(hkl, 2*(lmxh+1)**2)
        call dpzero(ghkl,3*(lmxh+1)**2)
      endif
C     lmax+1 required for gradients
C     lmax1 = lmax+1

C --- Start big loop over rsm-dependent l ---
C     Negative 1st smoothing radius => constant rsm
      fixdrs = rsm(0) < 0
      rsm0 = dabs(rsm(0))
      eh0 = eh(0)
      rsx = -9999
      ehx =  9999
C     This makes hanszd return true laplacian of hsm in ddhs
      jhd = 2 + 10*isw(ldot)
C     Decrement from highest l to simplify treatment of l-dependent rsm
      do  ilr = lmax, 0, -1
        if (.not. fixdrs) then
          rsm0 = dabs(rsm(ilr))
          eh0  = eh(ilr)
        endif
        if (dabs(rsm0-rsx) + dabs(eh0-ehx) > tol) then
          rsx = rsm0
          ehx = eh0

          call hanszd(jhd,r1,eh0,-rsm0,ilr+1,hs,dhs,ddhs,hsp,dhsp,ddhsp)
          if (ldot) then
            call dcopy(ilr+2,hsp,1,hs,1)
            call dcopy(ilr+2,ddhsp,1,ddhs,1)
          elseif (iopt == 3 .or. iopt == 4) then
            call radkj(eh0,r1,ilr+2,hsp,dhsp,dhsp,dhsp,0)
            call daxpy(ilr+2,-1d0,hsp,1,hs,1)
            call dscal(ilr+2,-eh0,hsp,1)
            call daxpy(ilr+2,-1d0,hsp,1,ddhs,1)
          endif
C         This would be more efficient, but gradient doesn't come out the same
C          if (iopt == 2 .or. iopt == 4) then
C            call radgkl(r1,-rsm0,1,ilr+2,1,gkl)
C            do  l = 0, ilr+1
C              hs(l) = hs(l) + cg0(l)*gkl(0,l)*r1**l
C              ddhs(l) = ddhs(l) + cg0(l)*gkl(1,l)*r1**l
C            enddo
C          endif

C     ... Solid functions and laplacians for given rsm0 up to curent lmax+1
          nlmp = (ilr+2)**2
          do  ilm = 1, nlmp
            l = ll(ilm)
            hloc(0,ilm) = hs(l)  *yl(ilm)
            hloc(1,ilm) = ddhs(l)*yl(ilm)
          enddo

C     ... Make grad hsm for l up to ilr = current lmax
          call dpzero(ghkl, 3*(ilr+1)**2)
          do  ilm = 1, (ilr+1)**2
            call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
            ghkl(1,ilm) = ghkl(1,ilm) -cx1*hloc(0,kx1)-cx2*hloc(0,kx2)
            ghkl(2,ilm) = ghkl(2,ilm) -cy1*hloc(0,ky1)-cy2*hloc(0,ky2)
            ghkl(3,ilm) = ghkl(3,ilm) - cz*hloc(0,kz)

            if (ilm <= ilr*ilr) then
              ghkl(1,kx1) = ghkl(1,kx1) - cx1*hloc(1,ilm)
              ghkl(1,kx2) = ghkl(1,kx2) - cx2*hloc(1,ilm)
              ghkl(2,ky1) = ghkl(2,ky1) - cy1*hloc(1,ilm)
              ghkl(2,ky2) = ghkl(2,ky2) - cy2*hloc(1,ilm)
              ghkl(3,kz)  = ghkl(3,kz)  - cz *hloc(1,ilm)
            endif
          enddo

c     ... Copy result to H_kL, up to current lmax only
          nlmp = (ilr+1)**2
          if (iopt <= 3) then
            call dcopy(2*(ilr+1)**2,hloc(0,1),1,hkl(0,1),1)
          endif

C     ... Add gaussians for functions that require it
          if (iopt == 2 .or. iopt == 4) then
            call solgsg(ph,-rsm0,ilr,1,1,gloc,ggrad)
            do  ilm = 1, nlmp
              l = ll(ilm)
              hkl(0,ilm) = hloc(0,ilm) + cg0(l)*gloc(0,ilm)
              hkl(1,ilm) = hloc(1,ilm) + cg0(l)*gloc(1,ilm)
              ghkl(1,ilm) = ghkl(1,ilm) + cg0(l)*ggrad(1,ilm)
              ghkl(2,ilm) = ghkl(2,ilm) + cg0(l)*ggrad(2,ilm)
              ghkl(3,ilm) = ghkl(3,ilm) + cg0(l)*ggrad(3,ilm)
            enddo
          endif

        endif
      enddo

      end
C
CC     Tests grad hsm.  grad j done through solhjg.
C      subroutine fmain
C      implicit none
C      integer lmxh,nlmh
C      parameter (lmxh=2, nlmh=(lmxh+1)**2)
C      double precision ph(3),php(3),eh(0:lmxh)
C      double precision hl(nlmh),jl(nlmh),ghl(3,nlmh),gjl(3,nlmh)
C      double precision hkl(0:1,nlmh),ghkl(3,nlmh),bl(nlmh)
C      double precision hlx(nlmh),blx(nlmh)
C      double precision sdiffh,sdiffj,xx,radius,ddot,gh,gj,rs
C      integer jlm,il,im
C
C      ph(1) = .3d0 *1
C      ph(2) = .9d0 *1
C      ph(3) = 1.3d0
C      eh = -.3d0 * 3
C      rs = 1.3d0
C
CC      ph(1) = 0.34212723867136507d0
CC      ph(2) = 0.24856998888378706d0
CC      ph(3) = 0.90617984593866396d0
CC      eh = -.7d0
C
CC      print *, '!!'
CC      radius = dsqrt(ddot(3,ph,1,ph,1))
CC      call dscal(3,1/radius,ph,1)
C
CC     Use solhjg for reference hl
C      call solhjg(ph,eh,lmxh,hlx,jl,ghl,gjl)
CC     Use soldhj for reference jl
C      call soldhj(ph,eh,0,lmxh,hl,bl)
CC     Testing this function:
C      call solhsg(ph,-.1d0,eh,lmxh,0,0d0,hkl,ghkl)
C
C      write (*,860) lmxh, eh(1), ph
C  860 format(/' Check H and J against other routines,',
C     .  '  lmxh =',i2,'  eh = ',f5.2,'  r =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  4x,'H,solhsg',4x,'H,solhjg',2x,
C     .  4x,'J,soldhj',4x,'J,solhjg',
C     .  4x,' diff(H)     diff(J)'/1x,72('-'))
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hkl(0,jlm)-hl(jlm))
C          sdiffj = sdiffj + dabs(bl(jlm)-jl(jlm))
C          write(*,870) jlm,il,im,hkl(0,jlm),hl(jlm),
C     .      bl(jlm),jl(jlm),
C     .      hkl(0,jlm)-hl(jlm),
C     .      bl(jlm)-jl(jlm)
C  870     format(i4,2i3,2f12.6,2x,2f12.6,1p,2e12.2)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C  871 format('  sum_L abs(diff) ',42x,1p,2e12.2)
C
C      write (*,861)
C  861 format(//' Check radial part of grad H and grad J against ',
C     .  'numerical differentiation'/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gH,num',3x,'gH,analytic',2x,
C     .  4x,'gJ,num',3x,'gJ,analytic',
C     .  2x,' diff(gH)    diff(gJ)'/1x,72('-'))
C
CC     Rewrite ph as radius * unit vector
C      radius = dsqrt(ddot(3,ph,1,ph,1))
C      call dscal(3,1/radius,ph,1)
C
C      php(1) =  (radius + .0001d0)*ph(1)
C      php(2) =  (radius + .0001d0)*ph(2)
C      php(3) =  (radius + .0001d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hlx,blx)
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      php(1) =  (radius - .0001d0)*ph(1)
C      php(2) =  (radius - .0001d0)*ph(2)
C      php(3) =  (radius - .0001d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
CC     hlx = (hlx - hl)/.0002d0
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx = (blx - bl)/.0002d0
C      call solhsg(radius*ph,-rs,eh,lmxh,0,0d0,hkl,ghkl)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do  im = -il, il
C          jlm = jlm + 1
C          gh = ddot(3,ghkl(1,jlm),1,ph,1)
C          gj = ddot(3,gjl(1,jlm),1,ph,1)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh,blx(jlm),gj,
C     .      hlx(jlm)-gh,blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C      write (*,862)
C  862 format(//' Check x component of grad H and grad J against ',
C     .  'numerical differentiation'/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gH,num',3x,'gH,analytic',2x,
C     .  4x,'gJ,num',3x,'gJ,analytic',
C     .  2x,' diff(gH)    diff(gJ)'/1x,72('-'))
C
C
C      php(1) =  (radius + .0000d0)*ph(1) + .0001d0
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hlx,blx)
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      php(1) =  (radius - .0000d0)*ph(1) - .0001d0
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx = (blx - bl)/.0002d0
C      call solhsg(radius*ph,-rs,eh,lmxh,0,0d0,hkl,ghkl)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do  im = -il, il
C          jlm = jlm + 1
C          gh = ghl(1,jlm)
C          gj = gjl(1,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh,blx(jlm),gj,
C     .      hlx(jlm)-gh,blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C      write (*,863)
C  863 format(//' Check y component of grad H and grad J against ',
C     .  'numerical differentiation'/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gH,num',3x,'gH,analytic',2x,
C     .  4x,'gJ,num',3x,'gJ,analytic',
C     .  2x,' diff(gH)    diff(gJ)'/1x,72('-'))
C
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2) + .0001d0
C      php(3) =  (radius + .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hlx,blx)
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2) - .0001d0
C      php(3) =  (radius - .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
CC     hlx = (hlx - hl)/.0002d0
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do  im = -il, il
C          jlm = jlm + 1
C          gh = ghl(2,jlm)
C          gj = gjl(2,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh,blx(jlm),gj,
C     .      hlx(jlm)-gh,blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C      write (*,864)
C  864 format(//' Check z component of grad H and grad J against ',
C     .  'numerical differentiation'/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gH,num',3x,'gH,analytic',2x,
C     .  4x,'gJ,num',3x,'gJ,analytic',
C     .  2x,' diff(gH)    diff(gJ)'/1x,72('-'))
C
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3) + .0001d0
C      call soldhj(php,eh,0,lmxh,hlx,blx)
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3) - .0001d0
C      call soldhj(php,eh,0,lmxh,hl,bl)
CC     hlx = (hlx - hl)/.0002d0
C      call solhsg(php,-rs,eh,lmxh,0,0d0,hkl,ghl)
C      hlx = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do  im = -il, il
C          jlm = jlm + 1
C          gh = ghl(3,jlm)
C          gj = gjl(3,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh,blx(jlm),gj,
C     .      hlx(jlm)-gh,blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC     Restore ph
C      call dscal(3,radius,ph,1)
C
C      end
