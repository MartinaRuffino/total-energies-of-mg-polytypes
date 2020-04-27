      subroutine sole0g(ph,rsm,eh,lmxh,iopt,kmax,c01,hkl,ghkl)
C- Solid smooth functions based sm Hankels, lapl and gradient, at point
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :coordinates relative to centre of function hkl
Ci         :NB: if 100s digit of iopt is nonzero, only ph(1) is used.
Ci   rsm   :vector of l-dependent smoothing radii of the function
Ci         :rsm<0 => rsm, eh are independent of l.  |rsm(0)| is used for all l
Ci         :rsm>0 => rsm,eh must be specified for 0..lmxh
Ci   eh    :l-dependent Hankel energies, but see rsm about l-dependence
Ci   lmxh  :Generate functions for l=0..lmxh
Ci   kmax  :0 => return function (hkl(0,:)) but not its Laplacian (hkl(1,:))
Ci         :1 => return function and its Laplacian both
Ci   iopt  :Determines the function to return; see Remarks
Ci         :1s digit fixes function.
Ci         :Note: kmax=0 => only function is returned
Ci         :      kmax=1 => Laplacian of function is also returned
Ci         :0 return h0l (and Laplacian if kmax=1, i,e. hkl for k=0:1)
Ci         :1 return energy derivative of h0l
Ci         :2 return h0l + c01(0)*g0l + c01(1)*g1l
Ci         :3 return h0l - h0l(rsm=0)
Ci         :4 return h0l - h0l(rsm=0) + c01(0)*g0l + c01(1)*g1l
Ci         :5 return g0l, g1l (NOT IMPLEMENTED)
Ci         :10s digit determines whether gradients are made
Ci         :0 Do not return gradients
Ci         :1 Return gradient of hkl(0:kmax,:)
Ci         :Gradients are returned in ghkl
Ci         :100s digit =>
Ci         :0 return solid functions and their derivatives.
Ci         :1 return only radial part of function and derivatives.
Ci         :  Program returns hkl(:,1:1+lmxh), not hkl(:,1:(1+lmxh)**2)
Ci         :  Also only radial derivative is made; ghkl is treated as
Ci         :  though it were dimensioned ghkl(0:kmax,1:1+lmxh).
Ci   c01   :Coefficients to amout of gkl to add to hkl; see Remarks
Co Outputs
Co   hkl   :makes function f at point ph for L = 1,...,(lmxh+1)^2
Co         :and, if kmax=1, nabla^2 f
Co         : mod(iopt,10)     function
Co         :  0               hsm(rsml,eh;ph)
Co         :  1               hsm-dot
Co         :  2               hsm + c01(0)*g0l + c01(1)*g1l
Co         :  3               hsm - hsm(rsm->0)
Co         :  4               hsm + c01(0)*g0l + c01(1)*g1l - hsm(rsm->0)
Co   ghkl  : mod(iopt,10)     function
Co         :  0               not touched
Co         :  1               \grad hkl, and grad Laplacian, if kmax=1
Cl Local variables
Cl  fixdrs :T => rsm are l-dependent
Cl  ghkll  : radial part of ghkl; used when 100s digit iopt is set
Cb Bugs
Cb   Laplacian of hdot is not correctly given (iopt=1,k=1)
Cb   ghkl is not returned when 100s digit iopt is set
Cr Remarks
Cr   Functions are designed for optimally chosen screened basis set.
Cr   Returns functions in hkl(0,:), Laplacians in hkl(1,:),
Cr   grad function in ghkl(1:3,0,:), grad Laplacian in ghkl(1;3,1,:)
Cr     Use iopt 0,1  for sm Hankel H0l or energy derivative
Cr              2-4  for hkL + c01(0)*g0l + c01(1)*g1l or similar.
Cr   Examples. Make:                                      iopt  kmax
Cr   h0l-h0l(rsm=0)                                         3     0
Cr   h0l+c01(0)*g0l+c01(1)*g1l, lap                         2     1
Cr   h0l+c01(0)*g0l+c01(1)*g1l-h0l(rsm=0), lap and grad    14     1
Cr
Cu Updates
Cu   19 Sep 11 Adapted from solhsg and solhjg
Cu   26 Aug 11 Additional kinds of functions may be calculated
Cu   20 Apr 11 Adapted from solgsg
C ----------------------------------------------------------------------
      implicit none
C Input parameters
      integer lmxh,iopt,kmax
      double precision ph(3),rsm(0:lmxh),eh(0:lmxh),c01(0:1,0:lmxh)
C Output parameters
      double precision hkl(0:kmax,(lmxh+1)**2),ghkl(3,0:kmax,*)
C Local variables
      logical fixdrs,ldot
      integer isw,l,ilm,nlm,ll,jhd,j,opt0,opt1,opt2
      double precision yl((lmxh+2)**2),gyl((lmxh+1)**2,3),px(3)
      double precision tol
      parameter (tol=1d-14)
      double precision rsm0,rsx,r1,ehx,eh0,ddot,dh,dlaph,cy1,xx,r1x
      double precision hs(0:lmxh),dhs(0:lmxh),ddhs(0:lmxh)
      double precision hsp(0:lmxh),dhsp(0:lmxh),ddhsp(0:lmxh)
      double precision dlh(0:lmxh)
      double precision ak(0:lmxh),aj(0:lmxh),dk(0:lmxh),dj(0:lmxh)
      double precision ddak(0:lmxh),dlak(0:lmxh)
      double precision gkl(0:2,0:lmxh),dgkl(0:2,0:lmxh),gkl2(0:2,0:lmxh)
      double precision ghkll(0:kmax,1:lmxh+1)

      opt0 = mod(iopt,10)
      opt1 = mod(iopt/10,10)
      opt2 = mod(iopt/100,10)
      ldot = opt0 == 1
      nlm  = (lmxh+1)**2
      call dpzero(dlh,1+lmxh)
      call dpzero(dlak,1+lmxh)

C --- Spherical harmonics and gradient (unit radius) ---
      if (opt2 /= 1) then
        r1 = dsqrt(ddot(3,ph,1,ph,1))
        if (r1 > 0) then
          call dpscop(ph,px,3,1,1,1/r1)
        else
          px(1) = 0
          px(2) = 0
          px(3) = 1
        endif
        call ropyln(1,px(1),px(2),px(3),lmxh+1,1,yl,xx)
        if (opt1 /= 0)
     .  call ropylg(1,lmxh,nlm,1,1,px(1),px(2),px(3),xx,yl,gyl)
        cy1 = dsqrt(3/(16*datan(1d0)))
      else
        r1 =  ph(1)
      endif

C --- Radial part of component functions and their derivatives ---
C     Negative 1st smoothing radius => constant rsm
      fixdrs = rsm(0) < 0
      rsm0 = dabs(rsm(0))
      eh0 = eh(0)
      rsx = -9999
      ehx =  9999
C     This makes hanszd return true laplacian of hsm in ddhs
      jhd = 2 + 10*isw(ldot)
C     Decrement from highest l to simplify treatment of l-dependent rsm
      do  l = lmxh, 0, -1
        if (.not. fixdrs) then
          rsm0 = dabs(rsm(l))
          eh0  = eh(l)
        endif
        if (dabs(rsm0-rsx) + dabs(eh0-ehx) > tol) then
C         Require derivative of Laplacian; do numerically for now
          if (kmax == 1 .and. opt1 == 1) then
            call hanszd(jhd,r1+1d-5,eh0,-rsm0,l,hs,dhs,dlh,
     .        hsp,dhsp,ddhsp)
            call hanszd(jhd,r1-1d-5,eh0,-rsm0,l,hs,dhs,ddhs,
     .        hsp,dhsp,ddhsp)
            dlh(0:l) = (dlh(0:l) - ddhs(0:l))/2d-5
          endif
          r1x = max(r1,1d-16)
          call hanszd(jhd,r1x,eh0,-rsm0,l,hs,dhs,ddhs,hsp,dhsp,ddhsp)
          if (ldot) then
            call dcopy(l+1,hsp,1,hs,1)
            call dcopy(l+1,dhsp,1,dhs,1)
            call dcopy(l+1,ddhsp,1,ddhs,1)
          endif
          if ((opt0 == 3 .or. opt0 == 4) .and.dabs(eh0-ehx) > tol) then
            call radkj(eh0,r1x,l,ak,aj,dk,dj,isw(ldot))
            call dpscop(ak,ddak,1+l,1,1,-eh0)
C           Require derivative of Laplacian
            if (kmax == 1 .and. opt1 == 1) then
              call dpscop(dk,dlak,1+l,1,1,-eh0)
            endif
          endif
          if (opt0 == 2 .or. opt0 == 4) then
            call radgkg(1,r1x,-rsm0,2,l,2,gkl,dgkl,gkl2)
          endif
          rsx = rsm0
          ehx = eh0
        endif
      enddo

C ... Assemble radial functions, derivatives from components
      if (opt0 >= 2) then
      do  l = 0, lmxh
        if (opt0 == 2 .or. opt0 == 4) then
          hs(l) =   hs(l) + c01(0,l)*gkl(0,l)  + c01(1,l)*gkl(1,l)
          dhs(l) =  dhs(l) + c01(0,l)*dgkl(0,l) + c01(1,l)*dgkl(1,l)
          ddhs(l) = ddhs(l) + c01(0,l)*gkl(1,l)  + c01(1,l)*gkl(2,l)
          dlh(l) =  dlh(l) + c01(0,l)*dgkl(1,l) + c01(1,l)*dgkl(2,l)
        endif
        if (opt0 == 3 .or. opt0 == 4) then ! subtract h(rs->0)
          hs(l) =   hs(l) - ak(l)
          dhs(l) =  dhs(l) - dk(l)
          ddhs(l) = ddhs(l) - ddak(l)
          dlh(l) =  dlh(l) - dlak(l)
        endif
      enddo
      endif

C --- Radial parts of functions and their derivatives ---
      if (opt2 /= 0) then
        do  l = 0, lmxh
          hkl(0,l+1) = hs(l)
          if (kmax == 1) hkl(1,l+1) = ddhs(l)
          if (opt1 /= 0) then
            ghkll(0,l+1) = dhs(l)
            if (kmax == 1) ghkll(1,l+1) = dlh(l)
          endif
        enddo
        if (opt1 /= 0) then
          call dcopy((kmax+1)*(lmxh+1),ghkll,1,ghkl,1)
        endif
        return
      endif

C --- Solid functions, and their gradients ---
      do  ilm = 1, nlm
        l = ll(ilm)
        hkl(0,ilm) = hs(l)*yl(ilm)
        if (kmax == 1) hkl(1,ilm) = ddhs(l)*yl(ilm)
        if (opt1 /= 0) then
C         Radial derivative, spherical coordinates
          dh = dhs(l)*yl(ilm)
          dlaph = dlh(l)*yl(ilm)
C         Split grad r- into x,y,z- components;
C         add contribution from grad YL
          do  j = 3, 1, -1
            xx = yl(j)/cy1                ! yl(2:3) are y,z
            if (j == 1) xx = yl(4)/cy1  ! yl(4) is x
            ghkl(j,0,ilm) = dh*xx + hs(l)*gyl(ilm,j)/r1
            if (kmax == 1) then
            ghkl(j,1,ilm) = dlaph*xx + ddhs(l)*gyl(ilm,j)/r1
            endif
          enddo
        endif
      enddo

      end
CC     Tests various branches of sole0g
C      subroutine fmain
C      implicit none
C      integer lmxh,nlmh
C      parameter (lmxh=2, nlmh=(lmxh+1)**2)
C      double precision ph(3),php(3),eh(0:lmxh),rsml(0:lmxh)
C      double precision c01(0:1,0:lmxh)
C      double precision hl(nlmh),hld(nlmh)
C      double precision hkl(0:1,nlmh),gkl(0:3,nlmh)
C      double precision dgkl(0:3,nlmh),ddgkl(0:3,nlmh)
C      double precision ghl(3,nlmh),ghld(3,nlmh)
C      double precision ghkl(3,0:1,nlmh),ghkld(3,0:1,nlmh)
C      double precision hlx(nlmh),dhlx(nlmh),blx(nlmh),dblx(nlmh),
C     .  hkld(0:1,nlmh)
C      double precision sdiffh,sdiffj,radius,ddot,gh,gj,xx
C      double precision hs(0:lmxh),dhs(0:lmxh),ddhs(0:lmxh)
C      double precision dddhs(0:lmxh),dgklp(0:3,nlmh)
CC     double precision hsp(0:lmxh),dhsp(0:lmxh),ddhsp(0:lmxh)
C      double precision ghkll(0:1,1:1+lmxh)
C      integer jlm,il,im
C
C      ph(1) = .3d0 *1
C      ph(2) = .9d0 *1
C      ph(3) = 1.3d0
C      eh = -.3d0 * 3
CC     rs = 1.3d0
C      rsml = (/1.5d0,1.4d0,1.3d0/)
C      c01(0,0:lmxh) = (/2.5702d0,2.1840d0,1.7309d0/)
C      c01(1,0:lmxh) = (/-0.2545d0,-0.1713d0,-0.1064d0/)
C
CC      ghkl = -99
CC      call sole0g(ph,rsml,eh,lmxh,10,1,c01,hkl,ghkl)
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
C
C
CC     Use solhjg for reference hl
C      call solhjg(ph,eh,lmxh,hlx,hld,ghl,ghld)
CC ... Testing this function:
C      call sole0g(ph,-.1d0,eh,lmxh,0,1,0d0,hkl,ghkl)
C      write (*,860) lmxh, eh(1), ph
C  860 format(/' Check H(rsm->0), lap against solhjg,',
C     .  '  lmxh =',i2,'  eh = ',f5.2,'  r =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  4x,'H,sole0g',4x,'H,solhjg',2x,
C     .  3x,'LH,sole0g',3x,'LH,solhjg',
C     .  4x,' diff(H)    diff(LH)'/1x,72('-'))
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hkl(0,jlm)-hlx(jlm))
C          sdiffj = sdiffj + dabs(hkl(1,jlm)+eh(il)*hlx(jlm))
C          write(*,870) jlm,il,im,hkl(0,jlm),hlx(jlm),
C     .      hkl(1,jlm),-eh(il)*hlx(jlm),
C     .      hkl(0,jlm)-hlx(jlm),
C     .      hkl(1,jlm)+eh(il)*hlx(jlm)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C  870 format(i4,2i3,2f12.6,2x,2f12.6,1p,2e12.2)
C  871 format('  sum_L abs(diff) ',42x,1p,2e12.2)
C
CC ... Testing this function:
C      write (*,760) lmxh, eh(1), ph, rsml
C  760 format(/' Check H, Hdot,',
C     .  '  lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  4x,'H,sole0g',4x,'H,hsmml',3x,
C     .  2x,'dot,sole0g',3x,'dot,hsmml',
C     .  3x,' diff(H)    diff(dot)'/1x,72('-'))
C      call hsmml(ph,rsml,eh,lmxh,hlx,blx)
C      call sole0g(ph,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      call sole0g(ph,rsml,eh,lmxh,1,1,0d0,hkld,ghkl)
C      sdiffh = 0d0
C      sdiffj = 0d0
C      jlm = 0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hkl(0,jlm)-hlx(jlm))
C          sdiffj = sdiffj + dabs(hkld(0,jlm)-blx(jlm))
C          write(*,870) jlm,il,im,
C     .      hkl(0,jlm),hlx(jlm),
C     .      hkld(0,jlm),blx(jlm),
C     .      hkl(0,jlm)-hlx(jlm),
C     .      hkld(0,jlm)-blx(jlm)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C
CC ... Testing this function:
C      write (*,660) lmxh, eh(1), ph, rsml
C  660 format(/' Check f=H+c0*g0+c1*g1, Lap,',
C     .  '  lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  4x,'f,sole0g',2x,'f,assembled',3x,
C     .  1x,'Lf,sole0g',2x,'Lf,assembled',
C     .  2x,' diff(f)    diff(Lf)'/1x,72('-'))
C      call solgsg(ph,rsml,lmxh,2,3,gkl,ghl)
C      call sole0g(ph,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      jlm = 0
C      do  il = 0, lmxh
C        do  im = -il, il
C          jlm = jlm+1
C          hlx(jlm) = hkl(0,jlm) +
C     .               c01(0,il)*gkl(0,jlm)+c01(1,il)*gkl(1,jlm)
C          blx(jlm) = hkl(1,jlm) +
C     .               c01(0,il)*gkl(1,jlm)+c01(1,il)*gkl(2,jlm)
C        enddo
C      enddo
C      call sole0g(ph,rsml,eh,lmxh,2,1,c01,hkl,ghkl)
C      sdiffh = 0d0
C      sdiffj = 0d0
C      jlm = 0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hkl(0,jlm)-hlx(jlm))
C          sdiffj = sdiffj + dabs(hkl(1,jlm)-blx(jlm))
C          write(*,870) jlm,il,im,
C     .      hkl(0,jlm),hlx(jlm),
C     .      hkl(1,jlm),blx(jlm),
C     .      hkl(0,jlm)-hlx(jlm),
C     .      hkl(1,jlm)-blx(jlm)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C      call info0(0,1,0,' Repeat, radial part only (iopt=112)')
C      radius = dsqrt(ddot(3,ph,1,ph,1))
CC      this one needs to be scaled by r^l
CC      do il = lmxh, 0, -1
CC        call radgkl(radius,rsml(il),2,il,3,gkl)
CC      enddo
C
CC      dgklp = 0; ddgkl = 0
CC      call radgkg(1,radius+1d-5,rsml,2,lmxh,3,gkl,dgklp,ddgkl)
CC      call radgkg(1,radius-1d-5,rsml,2,lmxh,3,gkl,dgkl,ddgkl)
CC      dgklp = (dgklp - dgkl)/2d-5
C
C      call hanszd(2,radius+1d-5,eh,rsml,lmxh,hs,dhs,dddhs,xx,xx,xx)
C      call hanszd(2,radius-1d-5,eh,rsml,lmxh,hs,dhs,ddhs,xx,xx,xx)
C      dddhs = (dddhs - ddhs)/2d-5
C      call radgkg(1,radius,rsml,2,lmxh,3,gkl,dgkl,ddgkl)
C      call hanszd(2,radius,eh,rsml,lmxh,hs,dhs,ddhs,xx,xx,xx)
C
CC      print *, '!!' ; c01 = 0
C
C      do  il = 0, lmxh
C        jlm = il+1
C        hlx(jlm) = hs(il) +
C     .             c01(0,il)*gkl(0,jlm)+c01(1,il)*gkl(1,jlm)
C       dhlx(jlm) = dhs(il) +
C     .             c01(0,il)*dgkl(0,jlm)+c01(1,il)*dgkl(1,jlm)
C        blx(jlm) = ddhs(il) +
C     .             c01(0,il)*gkl(1,jlm)+c01(1,il)*gkl(2,jlm)
C       dblx(jlm) = dddhs(il) +
C     .             c01(0,il)*dgkl(1,jlm)+c01(1,il)*dgkl(2,jlm)
C      enddo
C      call sole0g(radius,rsml,eh,lmxh,112,1,c01,hkl,ghkll)
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        jlm = il + 1
C        sdiffh = sdiffh + dabs(hkl(0,jlm)-hlx(jlm))
C        sdiffj = sdiffj + dabs(hkl(1,jlm)-blx(jlm))
C        write(*,872) il,
C     .    hkl(0,jlm),hlx(jlm),
C     .    hkl(1,jlm),blx(jlm),
C     .    hkl(0,jlm)-hlx(jlm),
C     .    hkl(1,jlm)-blx(jlm)
C      enddo
C      write(*,871) sdiffh,sdiffj
C  872 format(4x,i3,3x,2f12.6,2x,2f12.6,1p,2e12.2)
C
C      call info0(0,1,0,' Derivative, radial part only')
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        jlm = il + 1
C        sdiffh = sdiffh + dabs(ghkll(0,jlm)-dhlx(jlm))
C        sdiffj = sdiffj + dabs(ghkll(1,jlm)-dblx(jlm))
C        write(*,872) il,
C     .    ghkll(0,jlm),dhlx(jlm),
C     .    ghkll(1,jlm),dblx(jlm),
C     .    ghkll(0,jlm)-dhlx(jlm),
C     .    ghkll(1,jlm)-dblx(jlm)
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,560) lmxh, eh(1), ph, rsml
C  560 format(/' Check f1=H-H(rsm=0), f2=H+c0*g0+c1*g1-H(rsm=0),',
C     .  '  lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  4x,'f,sole0g',2x,'f,assembled',3x,
C     .  1x,'Lf,sole0g',2x,'Lf,assembled',
C     .  2x,' diff(f)    diff(Lf)'/1x,72('-'))
C
C      call sole0g(ph,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      call solhjg(ph,eh,lmxh,blx,hld,ghl,ghl)
C      call daxpy(nlmh,-1d0,blx,1,hkl,2)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      call sole0g(ph,rsml,eh,lmxh,2,1,c01,hkld,ghkl)
C      call daxpy(nlmh,-1d0,blx,1,hkld,2)
C      blx(1:nlmh) = hkld(0,1:nlmh)
C
C      call sole0g(ph,rsml,eh,lmxh,3,1,c01,hkl,ghkl)    ! f1
C      call sole0g(ph,rsml,eh,lmxh,4,1,c01,hkld,ghkl) ! f2
C
C      sdiffh = 0d0
C      sdiffj = 0d0
C      jlm = 0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hkl(0,jlm)-hlx(jlm))
C          sdiffj = sdiffj + dabs(hkld(0,jlm)-blx(jlm))
C          write(*,870) jlm,il,im,
C     .      hkl(0,jlm),hlx(jlm),
C     .      hkld(0,jlm),blx(jlm),
C     .      hkl(0,jlm)-hlx(jlm),
C     .      hkld(0,jlm)-blx(jlm)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C      write (*,"(/' Repeat with kmax=0')")
C      call sole0g(ph,rsml,eh,lmxh,3,0,c01,hl,ghkl)  ! f1
C      call sole0g(ph,rsml,eh,lmxh,4,0,c01,hld,ghkl) ! f2
C
C      sdiffh = 0d0
C      sdiffj = 0d0
C      jlm = 0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          sdiffh = sdiffh + dabs(hl(jlm)-hlx(jlm))
C          sdiffj = sdiffj + dabs(hld(jlm)-blx(jlm))
C          write(*,870) jlm,il,im,
C     .      hl(jlm),hlx(jlm),
C     .      hld(jlm),blx(jlm),
C     .      hl(jlm)-hlx(jlm),
C     .      hld(jlm)-blx(jlm)
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,861) 'radial part', lmxh, eh(1), ph, rsml
C  861 format(//' Check ',a,' of grad H and dot against ',
C     .  'numerical differentiation'/
C     .  ' lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gH,num',3x,'gH,analytic',2x,
C     .  2x,'gHdot,num',1x,'gHdot,analytic',
C     .  1x,'diff(gH)    diff(gJ)'/1x,72('-'))
C
CC     Rewrite ph as radius * unit vector
C      radius = dsqrt(ddot(3,ph,1,ph,1))
C      call dscal(3,1/radius,ph,1)
C
C      php(1) =  (radius + .0001d0)*ph(1)
C      php(2) =  (radius + .0001d0)*ph(2)
C      php(3) =  (radius + .0001d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = hkl(0,1:nlmh)
C      php(1) =  (radius - .0001d0)*ph(1)
C      php(2) =  (radius - .0001d0)*ph(2)
C      php(3) =  (radius - .0001d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,10,1,0d0,hkl,ghkl)
C      call sole0g(radius*ph,rsml,eh,lmxh,11,1,0d0,hkld,ghkld)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ddot(3,ghkl(1,0,jlm),1,ph,1)
C          gj = ddot(3,ghkld(1,0,jlm),1,ph,1)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,861) 'x component', lmxh, eh(1), ph, rsml
C
C      php(1) =  (radius + .0000d0)*ph(1) + .0001d0
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = hkl(0,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1) - .0001d0
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,10,1,0d0,hkl,ghkl)
C      call sole0g(radius*ph,rsml,eh,lmxh,11,1,0d0,hkld,ghkld)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(1,0,jlm)
C          gj = ghkld(1,0,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,861) 'y component', lmxh, eh(1), ph, rsml
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2) + .0001d0
C      php(3) =  (radius + .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = hkl(0,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2) - .0001d0
C      php(3) =  (radius - .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,10,1,0d0,hkl,ghkl)
C      call sole0g(radius*ph,rsml,eh,lmxh,11,1,0d0,hkld,ghkld)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(2,0,jlm)
C          gj = ghkld(2,0,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,861) 'z component', lmxh, eh(1), ph, rsml
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3) + .0001d0
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = hkl(0,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3) - .0001d0
C      call sole0g(php,rsml,eh,lmxh,0,1,0d0,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      call sole0g(php,rsml,eh,lmxh,1,1,0d0,hkl,ghkl)
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,10,1,0d0,hkl,ghkl)
C      call sole0g(radius*ph,rsml,eh,lmxh,11,1,0d0,hkld,ghkld)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(3,0,jlm)
C          gj = ghkld(3,0,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,862) 'x component', lmxh, eh(1), ph, rsml
C  862 format(/' Check ',a,' of f=Hsm',
C     .  ' and lap against numerical differentiation'/
C     .  ' lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gf,num',3x,'gf,analytic',2x,
C     .  2x,'gLf,num',3x,'gLf,analytic',
C     .  3x,'diff(gH)   diff(gLf)'/1x,72('-'))
C
C      php(1) =  (radius + .0000d0)*ph(1) + .0001d0
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      blx(1:nlmh) = hkl(1,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1) - .0001d0
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,0,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(1,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,10,1,c01,hkl,ghkl)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(1,0,jlm)
C          gj = ghkl(1,1,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,863) 'y component', lmxh, eh(1), ph, rsml
C  863 format(/' Check ',a,' of f1=H+c0*g0+c1*g1',
C     .  ' and lap against numerical differentiation'/
C     .  ' lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gf,num',3x,'gf,analytic',2x,
C     .  3x,'gLf,num',2x,'gLf,analytic',
C     .  3x,'diff(gH)   diff(gLf)'/1x,72('-'))
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2) + .0001d0
C      php(3) =  (radius + .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,2,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      blx(1:nlmh) = hkl(1,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2) - .0001d0
C      php(3) =  (radius - .0000d0)*ph(3)
C      call sole0g(php,rsml,eh,lmxh,2,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(1,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,12,1,c01,hkl,ghkl)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(2,0,jlm)
C          gj = ghkl(2,1,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC ... Testing this function:
C      write (*,864) 'z component', lmxh, eh(1), ph, rsml
C  864  format(/' Check ',a,' of f1=H+c0*g0+c1*g1-H(rsm=0)',
C     .  ' and lap against numerical differentiation'/
C     .  ' lmxh=',i1,'  eh=',f5.2,'  r =',3f7.4,'  rs =',3f7.4/
C     .  1x,72('-')/'   L  l  m',
C     .  5x,'gf,num',3x,'gf,analytic',2x,
C     .  3x,'gLf,num',2x,'gLf,analytic',
C     .  3x,'diff(gH)   diff(gLf)'/1x,72('-'))
C
C      php(1) =  (radius + .0000d0)*ph(1)
C      php(2) =  (radius + .0000d0)*ph(2)
C      php(3) =  (radius + .0000d0)*ph(3) + .0001d0
C      call sole0g(php,rsml,eh,lmxh,4,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = hkl(0,1:nlmh)
C      blx(1:nlmh) = hkl(1,1:nlmh)
C
C      php(1) =  (radius - .0000d0)*ph(1)
C      php(2) =  (radius - .0000d0)*ph(2)
C      php(3) =  (radius - .0000d0)*ph(3) - .0001d0
C      call sole0g(php,rsml,eh,lmxh,4,1,c01,hkl,ghkl)
C      hlx(1:nlmh) = (hlx(1:nlmh) - hkl(0,1:nlmh))/.0002d0
C      blx(1:nlmh) = (blx(1:nlmh) - hkl(1,1:nlmh))/.0002d0
C
C      call sole0g(radius*ph,rsml,eh,lmxh,14,1,c01,hkl,ghkl)
C
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghkl(3,0,jlm)
C          gj = ghkl(3,1,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
C      write (*,"(/' Repeat with kmax=0')")
C      call sole0g(radius*ph,rsml,eh,lmxh,14,0,c01,hl,ghl)
C      jlm = 0
C      sdiffh = 0d0
C      sdiffj = 0d0
C      do  il = 0, lmxh
C        do im = -il, il
C          jlm = jlm + 1
C          gh = ghl(3,jlm)
C          gj = ghkl(3,1,jlm)
C          sdiffh = sdiffh + dabs(hlx(jlm)-gh)
C          sdiffj = sdiffj + dabs(blx(jlm)-gj)
C          write(*,870) jlm,il,im,
C     .      hlx(jlm),gh, blx(jlm),gj,
C     .      hlx(jlm)-gh, blx(jlm)-gj
C        enddo
C      enddo
C      write(*,871) sdiffh,sdiffj
C
CC     Restore ph
C      call dscal(3,radius,ph,1)
C
C      end
