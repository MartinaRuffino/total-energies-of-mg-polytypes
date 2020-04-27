      subroutine ropylg(lp,lmax,ndim,nd,n,x,y,z,r2,yl,gyl)
C- Gradients of YL's (functions or polynomials) for a set of points, with YL as input
C ----------------------------------------------------------------------
Ci Inputs
Ci   lp    :0  => gradient of (r^l YL), YL = real harmonic function on (x,y,z)
Ci         :      gradient of Ylm should return linear combination of Y_l-1,m'
Ci         :1  => gradient of YL.
Ci         :      This consists of linear combinations of Y_l-1,m' and Y_l+1,m'
Ci         :-1 => gradient of (r^(-l-1) YL).
Ci         :      This consists of linear combinations of Y_l+1,m'
Ci   lmax  :maximum l
Ci   ndim  :dimensions gyl.  Must be at least (lmax+1)**2
Ci   nd    :leading dimension of yl,gyl
Ci   n     :number of points
Ci   x,y,z :cartesian coordinates of points
Ci   r2    :x^2+y^2+z^2
Ci   yl    :Spherical harmonic polynomials YL.  YL's must be normalized
Ci         :and generated through lmax+1 (i.e. nlm=1:(lmax+2)**2)
Co Outputs
Co   gyl   :gradient of yl or a gradient of (power of r) * yl
Cb Bugs
Cb   Modes not properly debugged for r2 different from one
Cr Remarks
Cr   r^p grad(r^-p YL) = grad(YL) + r^p grad(r^-p) YL
Cr   Second term is (-p*\mathbf{r}/r^2) YL
Cr
Cr   Note :  grad r^l*_Y_l yield combinations of  Y_l+1 only
Cr           grad r^-l-1*Y_l yield combinations of Y_l-1 only
Cr   See ylmbyr.
Cr
Cr   For debugging, see alternative to calculating grad YL using only CG coefficients.
Cu Updates
Cu   09 Apr 18  Worked out how grad can be calculated purely from CG coefficients
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lp,lmax,ndim,nd,n
      double precision x(n),y(n),z(n),yl(nd,*),gyl(nd,ndim,3),r2(n)
C ... Local parameters
      integer ilm,kx1,kx2,ky1,ky2,kz,l,m,i
      double precision cx1,cx2,cy1,cy2,cz,f1,f,fac
C ... debugging
C     integer nlm1
C     double precision coff(2,2,3,(lmax+1)**2)

      if ((lmax+1)**2 > ndim) call rx('ropylg: ndim too small')

C --- Gradients grad (r^l YL), or grad YL or grad (r^l-1 YL) ---
      ilm = 0
      do  l = 0, lmax
        fac = 2*l+1; if (lp>0) fac = l+1; if (lp<0) fac = 0
        do  m = -l, l
          ilm = ilm+1
          call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
          do  i = 1, n
            f1 = fac/r2(i)
            f  = (2*l+1)/r2(i)
            gyl(i,ilm,1) = yl(i,ilm)*x(i)*f1 - cx1*yl(i,kx1)*f - cx2*yl(i,kx2)*f
            gyl(i,ilm,2) = yl(i,ilm)*y(i)*f1 - cy1*yl(i,ky1)*f - cy2*yl(i,ky2)*f
            gyl(i,ilm,3) = yl(i,ilm)*z(i)*f1 - cz*yl(i,kz)*f
          enddo
        enddo
      enddo

C      call prmx('yl',yl,nd,n,ndim)
C      call prmx('gyl',gyl,nd,n,ndim*3)

C --- Debugging: alternative calculation of grad Y_L using only CG coefficients ---
C      if (lp /= 1) return
C
C      call dpzero(gyl,size(gyl))
C      nlm1 = lmax*lmax
C      ilm = 0
C      do  l = 0, lmax
C        do  m = -l, l
C          ilm = ilm+1
C          call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
CC         Y_l+1m' components of grad Y_l,m
C          f = l
C          do  i = 1, n
C            gyl(i,ilm,1) = gyl(i,ilm,1) - cx1*yl(i,kx1)*f - cx2*yl(i,kx2)*f
C            gyl(i,ilm,2) = gyl(i,ilm,2) - cy1*yl(i,ky1)*f - cy2*yl(i,ky2)*f
C            gyl(i,ilm,3) = gyl(i,ilm,3) - cz*yl(i,kz)*f
C          enddo
C
C          if (ilm > nlm1) cycle
C
CC         Y_lm' components of grad Y_l+1,m
C          f = l+2  ! f = ll(kx1)+1
C          do  i = 1, n
C            gyl(i,kx1,1) = gyl(i,kx1,1) + cx1*f*yl(i,ilm)
C            gyl(i,kx2,1) = gyl(i,kx2,1) + cx2*f*yl(i,ilm)
C            gyl(i,ky1,2) = gyl(i,ky1,2) + cy1*f*yl(i,ilm)
C            gyl(i,ky2,2) = gyl(i,ky2,2) + cy2*f*yl(i,ilm)
C            gyl(i,kz,3)  = gyl(i,kz,3)  + cz*f*yl(i,ilm)
C          enddo
C
C        enddo
C      enddo
CC     call prmx('gyl',gyl,nd,n,ndim*3)

      end

C#ifdefC TEST
CC     Test program to check ropylg
C      subroutine fmain
C      implicit none
C
C      integer mode,ilmx,nlm,nlmf,ip,lmax,j,ilm,l,lav,m,lp,pm
C      double precision pi,srfpi !,xx1,xx2
C      character*(120) strn,c*1
C      integer, parameter :: nnn=300
C      double precision p(3,nnn),wp(nnn)
C      integer np,nph,nth
C      procedure(integer) :: ll
C
C      real(8), allocatable :: xp(:),yp(:),zp(:),rp(:),r2(:)
C      real(8), allocatable :: yl(:,:),gyl(:,:,:),frl(:)
C      real(8), allocatable :: grp(:,:),ggrp(:,:),coff(:,:,:,:)
C
CC     mode=0  => show single YL on angular mesh can be decomposed into by integrating with Yl
CC     mode=1  => compare explicit polynomials to check normalization
CC     mode=2  => decompose (x,y,z)*Ylm into linear combination of Ylm
CC     mode=4  => decompose r^-l grad*(r^l Ylm) into linear combination of Ylm, Ylm = sph. harm. (no r^l)
CC     mode=14 => decompose      grad*(Ylm) into linear combination of Ylm
CC     mode=24 => decompose r^(l+1) grad*(r^(-l-1) Ylm) into linear combination of Ylm
CC     mode=8  => laplacians
C
C      mode = 14
CC     mode = 8
C      ilmx = 9
CC     ilmx = 2
C      lmax = 4
C      nlmf = (lmax+2)**2
C      nlm  = (lmax+1)**2
C      pi = 4*datan(1d0)
C      srfpi = dsqrt(4*pi)
C
CC --- Laplacian of Yl ---
C      nth = -122
C      call fpiint(nth,nph,np,p,wp)
C      call info2(1,0,0,' %i angular points.  mode=%i',np,mode)
C
C      allocate(xp(np),yp(np),zp(np),r2(np),rp(np),grp(np,3),ggrp(np,3))
C      allocate(frl(nlm),yl(np,nlmf),gyl(np,nlm,3))
C      allocate(coff(2,2,3,(lmax+1)**2))
C
CC     p(1:3,1) = [1d0,0d0,0d0]
C
C      call dcopy(np,p(1,1),3,xp,1)
C      call dcopy(np,p(2,1),3,yp,1)
C      call dcopy(np,p(3,1),3,zp,1)
C
C      call ropyln(np,xp,yp,zp,lmax+1,np,yl,r2)
CC     print *, '!! scale points on sphere by 10'
CC     xp = xp*10; yp = yp*10; zp = zp*10
C
C      if (mod(mode,10) == 0) then
C      call info2(1,1,0,' confirm integration of Yl(%i) on angular consists of Yl(%i)',ilmx,ilmx)
C      call xyl(nlm,np,wp,yl(1,ilmx),r2,yl,frl)
C      call pryl(0,' ... Yl proj Yl(r;l,m) for l=%i, m=%i: l,m,f = ',ilmx,1,frl,nlm,1d0)
C      endif
C
C      if (mod(mode,10) == 1) then
C      print *
C      print *, '... compare explicit polynomials to check normalization'
C      do  j = 1, 22, 5
C       print *, (3*zp(j)**2-1)*sqrt(5d0/4)/srfpi-yl(j,7) ! (3*z^2-1)*sqrt(5/4)/sqrt(4*pi)
C       print *, (xp(j)*zp(j))*sqrt(15d0/1)/srfpi-yl(j,8) ! (x^2-y^2)*sqrt(15/4)/sqrt(4*pi)
C       print *, (xp(j)**2-yp(j)**2)*sqrt(15d0/4)/srfpi-yl(j,9) ! (x^2-y^2)*sqrt(15/4)/sqrt(4*pi)
C       print *, yp(j)*(3*xp(j)**2-yp(j)**2)*sqrt(35d0/8)/srfpi-yl(j,10)
C       print *, xp(j)*yp(j)*zp(j)*sqrt(105d0)/srfpi-yl(j,11)
C       print *, yp(j)*(5*zp(j)**2-1)*sqrt(3*7d0/8)/srfpi-yl(j,12)
C       print *, zp(j)*(5*zp(j)**2-3)*sqrt(7d0/4)/srfpi-yl(j,13)
C       print *, xp(j)*(5*zp(j)**2-1)*sqrt(3*7d0/8)/srfpi-yl(j,14)
C       print *, zp(j)*(xp(j)**2-yp(j)**2)*sqrt(105d0/4)/srfpi-yl(j,15)
C       print *, xp(j)*(xp(j)**2-3*yp(j)**2)*sqrt(35d0/8)/srfpi-yl(j,16)
C      enddo
C      endif
C
C      if (mod(mode,10) == 2) then
C      print *
C      print *, '... decomposition of (x,y,z)*Ylm into linear combination of Ylm'
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = ilm-lav
C        if (j == 1) then
C          forall (j=1:np) rp(j) = yl(j,ilm)*xp(j)
C          call awrit2(' x*Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
C          pm = 1
C        endif
C        if (j == 2) then
C          forall (j=1:np) rp(j) = yl(j,ilm)*yp(j)
C          call awrit2(' y*Yl(%i%,2i): dl,d-m,cof= ',strn,len(strn),0,l,m)
C          pm = -1
C        endif
C        if (j == 3) then
C          forall (j=1:np) rp(j) = yl(j,ilm)*zp(j)
C          call awrit2(' z*Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
C          pm = 1
C        endif
C        call xyl(nlm,np,wp,rp,r2,yl,frl)
C        call pryl(0,strn,ilm,pm,frl,nlm,1d0)
C      enddo
C      enddo
C      endif
C
C      if (mod(mode,10) == 4) then
C      lp=0; if (mod(mode/10,10) == 1) lp=1; if (mod(mode/10,10) == 2) lp=-1
C      call ropylg(lp,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)  ! gradient and laplacian of spherical harmonic polynomials
C      call info0(1,1,0," ... CG coefficients C*Y_l'm' by tabulating grad*Ylm on angular and integrating (grad*Ylm) Yl'm'")
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = ilm-lav
C        if (m == -l) call info0(1,1,0,'')
C        call xyl(nlm,np,wp,gyl(1,ilm,j),r2,yl,frl)
C        if (j == 1) call awrit2(' gradx Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
C        if (j == 2) call awrit2(' grady Yl(%i%,2i): dl,d-m,cof= ',strn,len(strn),0,l,m)
C        if (j == 3) call awrit2(' gradz Yl(%i%,2i): dl,dm,cof = ',strn,len(strn),0,l,m)
C        pm = 1; if (j == 2) pm = -1
C        call pryl(0,strn,ilm,pm,frl,nlm,1d0)
C      enddo
C      enddo
C      endif
C
C      if (mod(mode,10) < 8) call rx0('done')
C
CC ... Show laplacian yl is 0:
C      call ropylg(0,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)  ! gradient and laplacian of spherical harmonic polynomials
C      call info2(1,1,0,
C     .  ' Gradient and Laplacian of spherical harmonic polynomials for ilm=%i',ilmx,ilmx)
C      call dpzero(ggrp,np*3)
C      do  j = 1, 3
C        call xyl(nlm,np,wp,gyl(1,ilmx,j),r2,yl,frl)
C        if (j == 1) c = 'x'
C        if (j == 2) c = 'y'
C        if (j == 3) c = 'z'
C        pm = 0
C        call pryl(0,' grad'//c//' Yl: l,m,cof = ',ilmx,pm,frl,nlm,1d0)
C        do  ilm = 1, nlm
C        do  ip = 1, np
C          ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
C        enddo
C        enddo
C        call xyl(nlm,np,wp,ggrp(1,j),r2,yl,frl)
C        call pryl(0,' '//c//' component of nabla',ilmx,1,frl,nlm,1d0)
C      enddo
C      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
C      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
C      call xyl(nlm,np,wp,ggrp,r2,yl,frl)
C      call pryl(0,' Laplacian: l,m,cof =',ilmx,0,frl,nlm,1d0)
CC     call cexit(1,1)
C
CCC ... Show laplacian r^-l yl is -l(l+1) yl;
CCC     Use grad(r^-1 yl) = 1/r sum_L (a_L yl)
CC      print *, 'gradient and Laplacian of spherical harmonics.  Should be -l(l+1) if r2=1'
CC      call ropylg(1,lmax,nlm,np,np,xp,yp,zp,r2,yl,gyl)
CCC ... Term 1/r grad r (grad r^-l yl)
CC      call dpzero(ggrp,np*3)
CC      do  20  j = 1, 3
CC        call xyl(nlm,np,wp,gyl(1,ilmx,j),r2,yl,frl)
CC        print *, j
CC        call prmr('component of grad',frl,nlm)
CC        do  18  ilm = 1, nlm
CC        do  18  ip = 1, np
CC   18   ggrp(ip,j) = ggrp(ip,j) + frl(ilm)*gyl(ip,ilm,j)
CC        call xyl(nlm,np,wp,ggrp(1,j),r2,yl,frl)
CCC        call prmr('1st term of nabla',frl,nlm)
CCC ... Term grad (1/r) . (grad r^-l yl)
CC        do  22  ip = 1, np
CC        if (j == 1) ggrp(ip,j) = ggrp(ip,j) - xp(ip)*gyl(ip,ilmx,j)
CC        if (j == 2) ggrp(ip,j) = ggrp(ip,j) - yp(ip)*gyl(ip,ilmx,j)
CC        if (j == 3) ggrp(ip,j) = ggrp(ip,j) - zp(ip)*gyl(ip,ilmx,j)
CC   22   continue
CC        call xyl(nlm,np,wp,ggrp(1,j),r2,yl,frl)
CCC        call prmr('component of nabla',frl,nlm)
CC   20 continue
CC      call dpadd(ggrp,ggrp(1,2),1,np,1d0)
CC      call dpadd(ggrp,ggrp(1,3),1,np,1d0)
CC      call xyl(nlm,np,wp,ggrp,r2,yl,frl)
CC      call prmr('laplacian',frl,nlm)
C
C      call cexit(1,1)
C      end
C
C      subroutine pryl(mode,strn,mlm,pm,fl,nlm,fac)
CC- Prints out nonzero (l,m) components of fl(ilm)*fac
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 print fac*fl
CCi         :1 print fac*sign(fl)*fl**2
CCi         :1 print fac*sign(fl)*fl**2*(2l+1)*l
CCi   strn
CCi   mlm   :ilm from which this function was derived
CCi   pm    :1   delta m = m(ilm) - m(fl)
CCi         :-1  delta m = m(ilm) - m(fl)
CCi         :0   print l and m not delta l and delta m
CCi   fl    :function as linear combation of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac
CCo Outputs
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   06 Apr 18
CC ----------------------------------------------------------------------
C      implicit none
C      character*(*) strn
C      integer mode,mlm,nlm,pm
C      double precision fl(nlm),fac
C      character*120 outs
C      integer ilm,l2,m2,lmax,l,m,lav,dl,dm
C      double precision f
C      real(8),parameter:: tol=1d-6
C      procedure(integer) :: ll
C
C      l = ll(mlm)
C      lav = l*l+l+1
C      m = mlm-lav
C      call awrit2(strn,outs,len(outs),0,l,m)
C
C      lmax = ll(nlm)
C      ilm = 0
C      do  l2 = 0, lmax
CC        ilm0 = l2**2 + l2 + 1   ! the m=0 element for this l
CC        ilmm = ilm0 - 2*l2      ! the m=0 element for this l-1
CC        ilmp = ilm0 + 2*(l2+1)  ! the m=0 element for this l+1
C        do  m2 = -l2, l2
C          ilm = ilm+1
C          if (abs(fl(ilm)) > tol) then
C            if (mode == 0) f = fl(ilm)*fac
C            if (mode == 1) f = fl(ilm)**2*fac*dsign(1d0,fl(ilm))
C            dl = l2-l; if (pm==0) dl = l2
C            dm = m2-m; if (pm<0) dm=m2+m; if (pm==0) dm = m2
C            call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,f)
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l-1))
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l+3))
C
CC           if (l2-l /= -1) stop 'oops'  ! True for grad r^l YL, YL = SH (not polynomials)
CC           if (l2-l /= 1) stop 'oops'  ! True for grad r^-l-1 YL, YL = SH (not polynomials)
C
C          endif
C        enddo
C      enddo
C      call info0(1,0,0,trim(outs))
C
C      end
C
C      subroutine xyl(nlm,np,wp,fp,r2,yl,fl)
CC- Yl-projection of function tabulated on an angular mesh
C      implicit none
C      integer nlm,np,ip,ilm
C      double precision fl(nlm),r2(np),fp(np),yl(np,nlm),wp(np)
C      double precision rl
C      procedure(integer) :: ll
C
C      call dpzero(fl,nlm)
C      do  ip = 1, np
C      do  ilm = 1, nlm
C        rl = dsqrt(r2(ip))**ll(ilm)
C        fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)/rl
C      enddo
C      enddo
C      end
C      subroutine prmr(strn,f,nl)
C      implicit none
C      integer nl,j,ifi
C      double precision f(nl)
C      character*(10) strn*(*), outs*80
C      ifi = 19
C      open(ifi,file='out')
C      write(ifi,334) nl, 2
C  334 format('% rows',i4,' cols', i2)
C      do  10  j = 1, nl
C   10 write(ifi,333) j, f(j)
C  333 format(i4, f18.10)
C      close(ifi)
C
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-6)
C      read(*,'(a80)') outs
C      if (outs == 'q') call rx0('quit in prmr')
C
C      end
C
C      subroutine makr(rsm,nr,x,y,z)
C      implicit none
C      integer nr,i,ir
C      double precision rs,rsm,x(1),y(1),z(1)
C      real ran1
C      rs = rsm
C      if (rsm < 1d-9) rs = .5d0
C      call ran1in(1)
C      do  10  i = 1, nr
C        ir = i+1
C        x(i) = abs((ran1()-.5d0)*5*rs)
C        y(i) = (ran1()-.5d0)*5*rs
C        z(i) = (ran1()-.5d0)*5*rs
C   10 continue
C
C      x(1) = .3d0*dsqrt(2d0)
C      y(1) = .4d0*dsqrt(2d0)
C      z(1) = .5d0*dsqrt(2d0)
C      end
C
C      subroutine prmx(strn,s,ns,nr,nc)
CC- writes matrix into out file (for debugging)
C      implicit none
C      integer nr,nc,ns,ifi
C      double precision s(ns,nc,2)
C      character*(14) fmt, fmt0, strn*(*), outs*80
C      integer i,j,fopna,i1mach
C      save fmt
C      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
C      ifi = fopna('out',29,0)
C      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#elseC
C      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
C#endifC
C      do  i = 1, nr
C        write (ifi,fmt) (s(i,j,1),j=1,nc)
C      enddo
C      write(ifi,*)
CC      do  12  i = 1, nr
CC   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
C      call fclose(ifi)
C
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#elseC
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endifC
C      read(*,'(a80)') outs
C
C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in prmx'
C#elseC
C      if (outs == 'q') call rx0('quit in prmx')
C#endifC
C      return
C
C      entry prmx0(fmt0)
C      fmt = fmt0
C      end
C
C#endif
