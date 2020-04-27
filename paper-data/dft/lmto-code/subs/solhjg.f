      subroutine solhjg(ph,eh,lmxh,hl,jl,ghl,gjl)
C- Real solid hankel and bessel functions, and gradients
C ----------------------------------------------------------------
Ci Inputs
Ci   ph    :(x,y,z) at which to evaluate hl,jl,ghl,gjl
Ci   eh    :Hankel energy
Ci   lmxh  :evaluate functions for L=1..(lmxh+1)**2
Co Outputs
Co   hl,jl: Hankel and Bessel functions for L=1..(lmxh+1)**2
Co   ghl  : gradient of hl
Co   gjl  : gradient of jl
Cr Remarks
Cr   MSM's standard defs, notes IV-43.
Cu Updates
Cu   12 Sep 11 Adapted from soldhj and gradfl
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer lmxh
      double precision ph(3),eh,hl(*),jl(*),ghl(3,*),gjl(3,*)
C Local parameters
      integer ilm,l,ll,nlm,j
      double precision xx,r1,ddot,cy1,dhloc,djloc,px(3)
      double precision yl((lmxh+2)**2),gyl((lmxh+1)**2,3)
      double precision ak(0:lmxh+1),aj(0:lmxh+1)
      double precision dk(0:lmxh+1),dj(0:lmxh+1)
      external ll

C     nlmp = (lmxh+2)**2
      nlm  = (lmxh+1)**2

C --- Spherical harmonics and gradient (unit radius) ---
      r1 = dsqrt(ddot(3,ph,1,ph,1))
      if (r1 > 0) call dpscop(ph,px,3,1,1,1/r1)
      call ropyln(1,px(1),px(2),px(3),lmxh+1,1,yl,xx)
      call ropylg(1,lmxh,nlm,1,1,px(1),px(2),px(3),xx,yl,gyl)
      cy1 = dsqrt(3/(16*datan(1d0)))

C     Radial part of gradient
      call radkj(eh,r1,lmxh,ak,aj,dk,dj,0)

C --- Solid Hankels, Bessels and their gradients ---
      do  ilm = 1, nlm
        l = ll(ilm)
        hl(ilm) = ak(l)*yl(ilm)
        jl(ilm) = aj(l)*yl(ilm)
C       Radial derivative, spherical coordinates
        dhloc = dk(l)*yl(ilm)
        djloc = dj(l)*yl(ilm)
C       Split grad r- into x,y,z- components;
C       add contribution from grad YL
        do  j = 3, 1, -1
          xx = yl(j)/cy1                ! yl(2:3) are y,z
          if (j == 1) xx = yl(4)/cy1  ! yl(4) is x
          ghl(j,ilm) = dhloc*xx + ak(l)*gyl(ilm,j)/r1
          gjl(j,ilm) = djloc*xx + aj(l)*gyl(ilm,j)/r1
        enddo
      enddo

      end

C      subroutine fmain
C      implicit none
C      integer lmxh,nlmh
C      parameter (lmxh=2, nlmh=(lmxh+1)**2)
C      double precision ph(3),php(3),eh(0:lmxh)
C      double precision hl(nlmh),jl(nlmh),ghl(3,nlmh),gjl(3,nlmh)
C      double precision hkl(0:1,nlmh),ghkl(3,nlmh),bl(nlmh)
C      double precision hlx(nlmh),blx(nlmh)
C      double precision sdiffh,sdiffj,radius,ddot,gh,gj
C      integer jlm,il,im
C
C      ph(1) = .3d0 *1
C      ph(2) = .9d0 *1
C      ph(3) = 1.3d0
C      eh = -.3d0 * 3
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
CC     Use solhsg for reference hl
C      call solhsg(ph,-.1d0,eh,lmxh,0,0d0,hkl,ghkl)
CC     Use soldhj for reference jl
C      call soldhj(ph,eh,0,lmxh,hl,bl)
CC     Testing this function:
C      call solhjg(ph,eh,lmxh,hlx,jl,ghl,gjl)
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
C
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
C      php(1) =  (radius - .0001d0)*ph(1)
C      php(2) =  (radius - .0001d0)*ph(2)
C      php(3) =  (radius - .0001d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
C      hlx = (hlx - hl)/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ddot(3,ghl(1,jlm),1,ph,1)
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
C      php(1) =  (radius - .0000d0)*ph(1) - .0001d0
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
C      hlx = (hlx - hl)/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
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
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2) - .0001d0
C      php(3) =  (radius - .0000d0)*ph(3)
C      call soldhj(php,eh,0,lmxh,hl,bl)
C      hlx = (hlx - hl)/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
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
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3) - .0001d0
C      call soldhj(php,eh,0,lmxh,hl,bl)
C      hlx = (hlx - hl)/.0002d0
C      blx = (blx - bl)/.0002d0
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
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
CC     call dscal(3,radius,ph,1)
C
C      end
