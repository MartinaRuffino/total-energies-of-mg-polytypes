      subroutine makus2(lpz,nr,rofi,g,gp,gz,phi,dphi,phip,dphip,phz,
     .  dphz,l,e,ez,z,v,ul,sl,ux,sx,ruu,rus,rss,ruz,rsz,rzz)
C- Kernel to make u and s from g and gp, for one l
C ----------------------------------------------------------------------
Ci Inputs
Cl   lpz   :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz
Cl         :1 value and slope of gz constructed to be zero at rmax
Cl         :  by admixture of phi,phidot
Cl         :2 a smooth Hankel tail is attached (extended local orbital)
Cl         :  but no alteration of gz is made
Cl         :3 same as lpz=2 for this routine.
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   g     :normalized wave function times r
Ci   gp    :energy derivative of g
Ci   gz    :(lpz=T) normalized semicore wave function times r
Ci   phi   :wave function at rmax, i.e. g/rmax
Ci   dphi  :slope of wave function at rmax, i.e. d(g/rmax)/dr
Ci   phip  :energy derivative of wave function at rmax
Ci   dphip :energy derivative of slope at rmax
Ci   phz   :value of sc w.f. at rmax: gz/rmax
Ci   dphz  :slope of sc wave function at rmax, i.e. d(gz/rmax)/dr
Ci   l     :l quantum number
Ci   e     :energy eigenvalue, needed for rel. small component
Ci   ez    :energy eigenvalue  of semicore state
Ci   z     :nuclear charge
Ci   v     :spherical potential
Co Outputs
Co   ul    :r * linear combination of wave functions; see Remarks
Co   sl    :r * linear combination of wave functions; see Remarks
Co   ux    :Analog of ul, but small component
Co   sx    :Analog of sl, but small component
Co   gz    :(lpz=T) returned as (input gz - gz(rmax) u - gz'(rmax) s)
Co   ruu   :diagonal product u_l*u_l,  including small component
Co   rus   :diagonal product u_l*s_l,  including small component
Co   rss   :diagonal product s_l*s_l,  including small component
Co   ...   The following are made when lpz=T:
Co   ruz   :diagonal product u_l*gz_l, including small component
Co   rsz   :diagonal product u_l*gz_l, including small component
Co   rzz   :diagonal product s_l*gz_l, including small component
Cr Remarks
Cr   This routine makes linear combinations (u,s) of out of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are as r * u and r * s, respectively.
Cu Updates
Cu   13 Jul 04 First implementation of extended local orbitals
Cu   21 Aug 01 extended to computation of semicore states
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: l,nr,lpz
      real(8), intent(in) :: dphi,dphip,e,ez,phi,phip,phz,dphz,z
      real(8), intent(in) :: g(nr*2),gp(nr*2),rofi(nr),v(nr)
      real(8), intent(inout) :: gz(nr*2), ul(nr), sl(nr), ruu(nr)
     . , rus(nr), rss(nr), ruz(nr), rsz(nr), rzz(nr), ux(nr), sx(nr)
C ... Local parameters
      integer :: ir,jr
      real(8) :: as,au,bs,bu,det,fllp1,gf11,gf12,gf22,r,tmcr,c,phzl,dphzl
      common /cc/ c

      fllp1 = l*(l+1)
      det = phi*dphip - dphi*phip
      au = dphip/det
      bu = -dphi/det
      as = -phip/det
      bs = phi/det
      if (lpz > 1) then
        phzl = 0
        dphzl = 0
      else
        phzl = phz
        dphzl = dphz
      endif

      if (z /= 0) then

C       This branch computes products of (g,gp,gz)
        do  ir = 1, nr
          jr = ir+nr
          r = rofi(ir)
          tmcr = r*c - (r*v(ir) - 2*z - r*e)/c
          gf11 = 1 + fllp1/tmcr**2

C         (u,s), large component

          ul(ir) = au*g(ir) + bu*gp(ir)
          sl(ir) = as*g(ir) + bs*gp(ir)
C         (u,s), small component
          ux(ir) = au*g(jr) + bu*gp(jr)
          sx(ir) = as*g(jr) + bs*gp(jr)

          ruu(ir) = gf11*ul(ir)*ul(ir) + ux(ir)*ux(ir)
          rus(ir) = gf11*ul(ir)*sl(ir) + ux(ir)*sx(ir)
          rss(ir) = gf11*sl(ir)*sl(ir) + sx(ir)*sx(ir)
        end do

        if (lpz /= 0) then
          do ir = 1, nr
            jr = ir+nr
            r = rofi(ir)
            tmcr = r*c - (r*v(ir) - 2*z - r*e)/c
            gf11 = 1 + fllp1/tmcr**2
            tmcr = r*c - (r*v(ir) - 2*z - r*ez)/c
            gf22 = 1 + fllp1/tmcr**2
            gf12 = (gf11 + gf22)*0.5d0

C           Subtract (phz ul + dphz sl) from gz
            gz(ir) = gz(ir) - phzl*ul(ir) - dphzl*sl(ir)
            gz(jr) = gz(jr) - phzl*ux(ir) - dphzl*sx(ir)

            ruz(ir) = gf12*ul(ir)*gz(ir) + ux(ir)*gz(jr)
            rsz(ir) = gf12*sl(ir)*gz(ir) + sx(ir)*gz(jr)
            rzz(ir) = gf22*gz(ir)*gz(ir) + gz(jr)*gz(jr)
          enddo
        end if

C --- Treat z=0 nonrelativistically ---
      else
        do  ir = 1, nr
          r = rofi(ir)
          ul(ir) = au*g(ir) + bu*gp(ir)
          sl(ir) = as*g(ir) + bs*gp(ir)

          ruu(ir) = ul(ir)*ul(ir)
          rus(ir) = ul(ir)*sl(ir)
          rss(ir) = sl(ir)*sl(ir)

          if (lpz == 0) cycle

C         Subtract (phzl ul + dphzl sl) from gz
          gz(ir) = gz(ir) - phzl*ul(ir) - dphzl*sl(ir)

          ruz(ir) = ul(ir)*gz(ir)
          rsz(ir) = sl(ir)*gz(ir)
          rzz(ir) = gz(ir)*gz(ir)
        enddo
      endif

C      call prrmsh('ruu',rofi,ruu,nr,nr,1)
C      call prrmsh('rus',rofi,rus,nr,nr,1)
C      call prrmsh('rss',rofi,rss,nr,nr,1)
C      if (lpz) then
C        call prrmsh('ul',rofi,ul,nr,nr,1)
C        call prrmsh('sl',rofi,sl,nr,nr,1)
C        call prrmsh('gz',rofi,gz,nr,nr,2)
C         call prrmsh('ruz',rofi,ruz,nr,nr,1)
C         call prrmsh('rsz',rofi,rsz,nr,nr,1)
C         call prrmsh('rzz',rofi,rzz,nr,nr,1)
C      endif
      end

      subroutine makus3(lpz,nr,rofi,l,e1,ez1,e2,ez2,z,v,
     .  ul1,sl1,ux1,sx1,gz1,ul2,sl2,ux2,sx2,gz2,
     .  ruu,rus,rss,ruz,rsz,rzz,rsu,rzu,rzs)
C- uu,us,ss products for crossed spins
C ----------------------------------------------------------------------
Ci Inputs
Cl   lpz   :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz
Cl         :1 value and slope of gz constructed to be zero at rmax
Cl         :  by admixture of phi,phidot
Cl         :2 a smooth Hankel tail is attached (extended local orbital)
Cl         :  but no alteration of gz is made
Cl         :3 same as lpz=2 for this routine.
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   g     :normalized wave function times r
Ci   gp    :energy derivative of g
Ci   gz    :(lpz=T) normalized semicore wave function times r
Ci   phi   :wave function at rmax, i.e. g/rmax
Ci   dphi  :slope of wave function at rmax, i.e. d(g/rmax)/dr
Ci   phip  :energy derivative of wave function at rmax
Ci   dphip :energy derivative of slope at rmax
Ci   phz   :value of sc w.f. at rmax: gz/rmax
Ci   dphz  :slope of sc wave function at rmax, i.e. d(gz/rmax)/dr
Ci   l     :l quantum number
Ci   e     :energy eigenvalue, needed for rel. small component
Ci   ez    :energy eigenvalue  of semicore state
Ci   z     :nuclear charge
Ci   v     :spherical potential
Co Outputs
Co   ul    :r * linear combination of wave functions; see Remarks
Co   sl    :r * linear combination of wave functions; see Remarks
Co   gz    :(lpz=T) returned as (input gz - gz(rmax) u - gz'(rmax) s)
Co   ruu   :diagonal product u_l*u_l,  including small component
Co   rus   :diagonal product u_l*s_l,  including small component
Co   rss   :diagonal product s_l*s_l,  including small component
Co   ...   The following are made when lpz=T:
Co   ruz   :diagonal product u_l*gz_l, including small component
Co   rsz   :diagonal product u_l*gz_l, including small component
Co   rzz   :diagonal product s_l*gz_l, including small component
Cr Remarks
Cr   This routine makes linear combinations (u,s) of out of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are as r * u and r * s, respectively.
Cu Updates
Cu   19 Dec 13 Adapted from makus2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: l,nr,lpz
      real(8) :: e1,ez1,e2,ez2,z
      real(8), intent(in) :: rofi(nr),v(nr),
     .  ul1(nr),sl1(nr),gz1(nr*2),ux1(nr),sx1(nr),
     .  ul2(nr),sl2(nr),gz2(nr*2),ux2(nr),sx2(nr)
      real(8), intent(inout) ::
     .  ruu(nr),rus(nr),rss(nr),ruz(nr),rsz(nr),rzz(nr),
     .  rsu(nr),rzu(nr),rzs(nr)
C ... Local parameters
      integer ir,jr
      real(8) :: fllp1,gf11,gf12,gf22,r,tmcr,c
      common /cc/ c

      fllp1 = l*(l+1)

      if (z /= 0) then

        do  ir = 1, nr
          jr = ir+nr
          r = rofi(ir)
          tmcr = r*c - (r*v(ir) - 2d0*z - r*(e1+e2)/2)/c
          gf11 = 1d0 + fllp1/tmcr**2

          ruu(ir) = gf11*ul1(ir)*ul2(ir) + ux1(ir)*ux2(ir)
          rus(ir) = gf11*ul1(ir)*sl2(ir) + ux1(ir)*sx2(ir)
          rsu(ir) = gf11*sl1(ir)*ul2(ir) + sx1(ir)*ux2(ir)
          rss(ir) = gf11*sl1(ir)*sl2(ir) + sx1(ir)*sx2(ir)

          if (lpz == 0) cycle

          tmcr = r*c - (r*v(ir) - 2d0*z - r*(ez1+ez2)/2)/c
          gf22 = 1d0 + fllp1/tmcr**2
          gf12 = (gf11 + gf22)/2

          ruz(ir) = gf12*ul1(ir)*gz2(ir) + ux1(ir)*gz2(jr)
          rzu(ir) = gf12*gz1(ir)*ul2(ir) + gz1(jr)*ux2(ir)
          rsz(ir) = gf12*sl1(ir)*gz2(ir) + sx1(ir)*gz2(jr)
          rzs(ir) = gf12*gz1(ir)*sl2(ir) + gz1(jr)*sx2(ir)
          rzz(ir) = gf22*gz1(ir)*gz2(ir) + gz1(jr)*gz2(jr)
        enddo

C --- Treat z=0 nonrelativistically ---
      else
        do  ir = 1, nr
          r = rofi(ir)

          ruu(ir) = ul1(ir)*ul2(ir)
          rus(ir) = ul1(ir)*sl2(ir)
          rsu(ir) = sl1(ir)*ul2(ir)
          rss(ir) = sl1(ir)*sl2(ir)

          if (lpz == 0) cycle

          ruz(ir) = ul1(ir)*gz2(ir)
          rzu(ir) = gz1(ir)*ul2(ir)
          rsz(ir) = sl1(ir)*gz2(ir)
          rzs(ir) = gz1(ir)*sl2(ir)
          rzz(ir) = gz1(ir)*gz2(ir)
        enddo
      endif

C      call prrmsh('ruu',rofi,ruu,nr,nr,1)
C      call prrmsh('rus',rofi,rus,nr,nr,1)
C      call prrmsh('rss',rofi,rss,nr,nr,1)
C      if (lpz) then
C        call prrmsh('ul',rofi,ul,nr,nr,1)
C        call prrmsh('sl',rofi,sl,nr,nr,1)
C        call prrmsh('gz',rofi,gz,nr,nr,2)
C         call prrmsh('ruz',rofi,ruz,nr,nr,1)
C         call prrmsh('rsz',rofi,rsz,nr,nr,1)
C         call prrmsh('rzz',rofi,rzz,nr,nr,1)
C      endif
      end
