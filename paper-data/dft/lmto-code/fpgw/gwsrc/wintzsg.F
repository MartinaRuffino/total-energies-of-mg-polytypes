      complex(8) function wintzsg_npm(npm,v,v0,x,wt,a,expa,we,
     i                                nx,esmr)
C- gaussian smeared-pole version of wintz.
C takao developed from wintz by ferdi.
C Assume that each eigenvale, and w-e as a result, has the width of esmr( Ry).
c
      implicit none
      integer(4)::npm
      integer(4) :: nx,i
      real(8) :: x(nx),wt(nx),expa(nx),a,we,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss,sig2,omd2w
      complex(8) ::v(npm*nx),v0,sum,   wintzsg1  !, wintzsg2,sum2
      complex(8) :: wwz,wintz_npm, img=(0d0,1d0)
      real(8)    :: esmr,wintsf,errsum=0d0,derfcx,derfc
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0
      !rmax =2 is by ferdi. Is it optimum? See wintz
      real(8)::sig, smxowx, ee, omg, ww,cons1,cons2,xx,aw,eee,aw2
      integer(4) ix,verbose
      logical :: test
c      real(8)::sumgauss,sumgauss1,sumgauss2 ! I found this bug 9Aug2004
      complex(8)::sumgauss,sumgauss1,sumgauss2
c      integer(4):: version

      logical ::  test_symmetric_W

c---------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' wintzsg: start!'
c      sig= 1d0 ! sig is in Ry.
c      ee = 2d0
c      do ix=1,2
c        omg= 0.1d0*(ix-1)
c        write(6,"(' omg smear=',f8.3, 2d23.15)")
c     &  omg, smxowx(omg,ee,sig), ee/(omg**2+ee**2)
c! smeared ee/(ww^2 + ee^2)
c      enddo
c      stop ' wintzsg: test end----------'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      we=1d0
c      do  i = 1,50
c        xx   = 1d-6*2**i
c        print *, ' i=',i,xx
c        sig2 =0.01d0
c        omd2w= sig2 *xx
c        cons1=(1d0 - exp (- omd2w/sig2))*we/omd2w
c        cons2 = we*( 1d0/sig2 - omd2w/sig2**2/2d0
c     &               + omd2w**2/sig2**3/6d0 - omd2w**3/sig2**4/24d0
c     &               + omd2w**4/sig2**5/120d0
c     &               - omd2w**5/sig2**6/720d0 )
c        write(6,"(' cons1=',3d25.16)")xx,cons1
c        write(6,"(' cons2=',3d25.16)")xx,cons2, cons2-cons1
c        write(6,*)
c      enddo
c      stop 'test end'


c- check part by direct integration by wintz.
c- this procedure is equivalent to whole of this routine. But time-consuming.
c      wwz =0d0
c      do ie = 1, nav
c        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
c        www=1d0
c        if(ie==1.or.ie==nav) www=0.5d0
c        wwz = wwz+ www*wintz(v,v0,x,wt,a,expa,wex,nx)
cxx       write(1112,"(4d24.16)") wex,wintz(v,v0,x,wt,a,expa,wex,nx)
c      enddo
c      wintzsg1=wwz/(nav-1)
c      wintzsg=wintzsg1
c      return

      if(esmr==0d0) then
        wintzsg_npm = wintz_npm(npm,v,v0,x,wt,a,expa,we,nx)
        return
      endif
      sig  = .5d0*esmr
      sig2 = 2d0*sig**2
      we2  = we**2
      sum  = 0d0

C --- simple integration scheme.
      do  i = 1,nx
        omd   = 1d0/x(i) - 1d0
        omd2w = omd**2 + we2
! pole weight is given by   1- exp (-R^2/sig2) = 1/N \int_0^R exp(- x^2/sig2) 2 pi r dr
! Gauss theorem---> but is this correct? Not three dimentional...
        if(omd2w/sig2 > 5d-3) then
          cons = (1d0 - exp (- omd2w/sig2))/omd2w
        else
          cons = ( 1d0/sig2 - omd2w/sig2**2/2d0
     &               + omd2w**2/sig2**3/6d0  - omd2w**3/sig2**4/24d0
     &               + omd2w**4/sig2**5/120d0- omd2w**5/sig2**6/720d0 )
        endif
        sum  = sum + we*cons*( v(i)-v0*expa(i) ) *wt(i)/(x(i)**2)
        if(test_symmetric_W()) then
        else
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)/(x(i)**2)
          endif
        endif
      enddo

C --- Gaussian part -------------------------------------
c      if(version()==-1) then !       print *,' version=',version()
c        sumgauss=
c     &          - .5d0*v0*dexp(we2*a**2)
c     &            *( derfcx(a*we) - derfcx(sqrt(a**2+1d0/sig2)*we) )
c      else
        aw = abs(a*we)
        aw2 = aw**2
        if(aw<6d0) then
          eee = we**2/sig2
          sumgauss = dsign(1d0,we)*.5d0*v0*exp(aw2)*( derfcx(sqrt(aw2+eee)) -derfcx(aw) )
c        sumgauss1=sumgauss
        else !April2004
          sumgauss =  -.5d0*v0* (2d0/(a*we)/pi)*
     &    sqrt(pi)/2d0*
     &       (1d0 - 1d0/2d0/aw2*
     &         (1d0 - 3d0/2d0/aw2*
     &           (1d0 - 5d0/2d0/aw2*
     &             (1d0 - 7d0/2d0/aw2*
     &               (1d0 - 9d0/2d0/aw2*
     &                 (1d0 - 11d0/2d0/aw2
     &     ))))))
cc       sumgauss2=sumgauss
c        endif
      endif
      wintzsg_npm = -sum/pi  + sumgauss
      if(verbose()>90)
     & write(1116, "(' we sig wintzsg_npm =',2f8.3,4f14.6)" )we,sig,wintzsg_npm
      return
c---test code
      print *,' sum= ', sum
      print *,' awe1= ', a*we
      print *,' erf1= ', derfcx(a*we)
      print *,' awe2= ', sqrt(a**2+1d0/sig2)*we
      print *,' erf2= ', derfcx(sqrt(a**2+1d0/sig2)*we)
      print *,' dexp= ', dexp(we2*a**2)
      print *,' erf3= ', derfcx(aw)
      print *,' erf4= ', derfcx(sqrt(aw2+eee))
      print *,' aw2= ', aw2
      aw = a*we
      write(6, "(' wintzsg_npm chk1=',2f8.3, 4d18.10)" )we,sig,
     &  - .5d0*v0*  dexp(aw**2)*( derfcx(aw) - derfcx(sqrt(a**2+1d0/sig2)*we)  )
      write(6, "(' wintzsg_npm chk2=',2f8.3, 4d18.10)" )we,sig,sumgauss1
      write(6, "(' wintzsg_npm chk3=',2f8.3, 4d18.10)" )we,sig,sumgauss2
      return
      end




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------------
      real(8) function smxowx(ww,ee,sig)
      implicit none
      real(8),parameter:: pi=3.1415926535897932d0
      real(8):: sq2pi,s2int
      real(8):: sum, fac,wwsig,s2,ww,ee,sig,dn2fac
      integer(4)::  mxn,n
      sq2pi=sqrt(2d0*pi)
      fac= exp( - ee**2/sig**2/2d0)*ee/sig
      mxn= 30 !???
      wwsig= ww/sig
      s2  = s2int(wwsig)
      print *,'wwsig=',wwsig
      print *,' s2=',s2
      sum = fac*s2
      dn2fac= 1d0
      do n=2,mxn
        fac = fac*(ee/sig)**2 /(2*n-1)/(2*n-2)
        dn2fac = dn2fac*(2*n-1)
        s2  = sq2pi*dn2fac - wwsig**2*s2
        sum = sum + fac*s2
        print *,' n sum',n,sum/sqrt(2*pi)
      enddo
      smxowx=sum/sqrt(2d0*pi)/sig
      end
      real(8) function s2int(beta)
      implicit none
      real(8),parameter:: pi=3.1415926535897932d0
      real(8):: derfc,sqpi,beta,derf
      real(8):: mu
      sqpi= sqrt(pi)
      mu  = 1d0/sqrt(2d0)
c      print *, '----------------beta*mu=',beta*mu
      s2int = sqpi/mu - pi*beta*exp(beta**2/2d0)*derfc(beta*mu)
      end


c -------------------------------------------------------------------
c --------------------------------------------------------------------
      complex(8) function wintzav(v,v0,x,wt,a,expa,we,
     i                                nx,esmr)
C- Smeared-pole version of wintz. takao developed from wintz by ferdi.
C Assume thaat each eigenvale, and w-e as a result, has the width of esmr( Ry).
c
      implicit none
      integer(4) :: nx,i
      real(8) :: x(nx),wt(nx),expa(nx),a,we,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss
      complex(8) ::v(nx),v0,sum,   wintzav1  !, wintzav2,sum2
      complex(8) :: wwz,wintz
      real(8)    :: esmr,wintsf,errsum=0d0
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
     &
c---------------------------------------------------
c      print *,' wintzav:'
      if(esmr==0d0) then
        wintzav=wintz(v,v0,x,wt,a,expa,we,nx)
        return
      endif
c-----------------------------------------------------
c- check part by direct integration by wintz.
c- this procedure is equivalent to whole of this routine. But time-consuming.
c      wwz =0d0
c      do ie = 1, nav
c        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
c        www=1d0
c        if(ie==1.or.ie==nav) www=0.5d0
c        wwz = wwz+ www*wintz(v,v0,x,wt,a,expa,wex,nx)
cxx       write(1112,"(4d24.16)") wex,wintz(v,v0,x,wt,a,expa,wex,nx)
c      enddo
c      wintzav1=wwz/(nav-1)
c      wintzav=wintzav1
c      return

c frequency integral
      we2  = we**2
      weh  = we + 0.25d0*esmr !because we is in a.u and esmr is in Ry.
      wel  = we - 0.25d0*esmr
      weh2 = weh**2
      wel2 = wel**2
      sum  = 0d0
c
c      sum2 = 0d0
c      wintzav2 = 0d0

ccase1 --------------------
      if ( max(dabs(weh),dabs(wel)) < rmax/a) then
ccase2 --------------------
ccase2 if (.false.) then
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
c   x = 1/(1+w')
c   We smear each pole with a pole distribution between wel ans weh,
c   So  it means \int_wel^weh dwe [ wt(i)/(x(i)**2) * we/(we2 + omd**2)]
          cons  = 0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum   = sum + wt(i)*(v(i) - v0*expa(i))*cons
ccccccccccccccccccccccccccccccccccccccccccccccccc
c        sum2  = sum2 + wt(i)*v(i) *cons
c      write(6,"(d15.6,2x,2d15.6,2x,2d15.6)")
c     & x(i), v(i)*cons,v0*expa(i)*cons
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & x(i), v(i),v0*expa(i)
ccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo

        if(abs(weh)>abs(wel)) then
          sss = 1d0
          rup = abs(weh)
          rdn = abs(wel) !;    iix=1
        else
          sss = -1d0
          rup = abs(wel)
          rdn = abs(weh) !;   iix=11
        endif
        wintzav  = - sum/pi - 0.5d0*v0* sss* wintsf(a*rdn,a*rup)/a
c       wintzav2 = - sum2/pi
ccccccccccccccccccccccccccccccccccccccccccccccccc
c      xxx=0d0
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & xxx, v0 ,v0
c      write(6,"(' - sum/pi =',2d13.6)") -sum/pi
c      write(6,"(' wintzav  =',2d13.6)")   wintzav
c      write(6,"(' - sum2/pi =',2d13.6)") -sum2/pi
c      stop ' ==== test end of wintz ========'
ccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c simple integration scheme.
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
          cons =  0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum  = sum + wt(i)*v(i)*cons
        enddo
        wintzav = -sum/pi  !;   iix =2
      endif
      wintzav = wintzav/(0.5d0*esmr) !because esmr is in Ry.

cccccccccccccccccccccccccccc
c      write(1116, "(' wel weh wintzav =',2f8.3,4f14.6)" )wel,weh,wintzav
ccccccccccccccccccccccccccc

c      wintzav2 = wintzav2/(0.5d0*esmr)
c      write(1116, "(' wintzav diff =',4d13.6,f10.6)" ) wintzav,wintzav2
c     &   ,abs(wintzav2-wintzav)/abs(wintzav)
cccccccc
c      stop ' ==== test end of wintz ========'
c      if(wintzav/=0d0 ) then
c        err = abs(wintzav-wintzav1)
c        if(err>1d-6) then
c          write(3302,"(3d24.16, i3,' err= ',2d13.5)")
c     &    we, wintzav, iix, err
c        endif
c     endif
cccccccccccccccccccccccccccc
      end

c----------------------------------------------------
      real(8) function wintsf(rdn,rup)
c integral \int^rup_rdn exp(x**2)*derfc(x)
c up to x <rmax<15 which is limited by the range of derfc(x).
      implicit none
      integer:: ix
      integer,parameter:: nd=2001
      real(8) :: polinta, rmax=14.9d0,rdn,rup,sumx,x,x0,xm,derfc
      real(8),save :: eeint(nd),xval(nd)
      logical :: init=.true.
c--------------
Cstop2rx 2013.08.09 kino      if(rup>rmax) stop 'wintsf: too large rup'
      if(rup>rmax) call rx( 'wintsf: too large rup')
Cstop2rx 2013.08.09 kino      if(rup<rdn ) stop 'wintsf: rup<rdn'
      if(rup<rdn ) call rx( 'wintsf: rup<rdn')
      if(init) then
        eeint(1) = 0d0
        xval(1)  = 0d0
        do ix= 2,nd
          x0= rmax*dble(ix-2)/(nd-1)
          x = rmax*dble(ix-1)/(nd-1)
          xm= .5d0*(x + x0)
          xval(ix) = x
          sumx  = 1/6d0 * exp(x0**2)*derfc(x0)
     &         + 4/6d0 * exp(xm**2)*derfc(xm)
     &         + 1/6d0 * exp(x **2)*derfc(x )
          eeint(ix) = eeint(ix-1) + rmax/(nd-1) * sumx
ccccccccccccccccccccccc
c        write(1113,"(4d24.16)") x,derfc(x),exp(x**2)*derfc(x),eeint(ix)
ccccccccccccccccccccccc
        enddo
        init=.false.
      endif
      wintsf = polinta(rup,xval,eeint,nd) - polinta(rdn,xval,eeint,nd)
      end


c--------------------------------------------------------------------
c--------------------------------------------------------------------
      complex(8) function wintz_npm(npm,v,v0,x,wt,a,expa,we,
     d                                nx)
c takao complex version of wint by ferdi
c 92.02.24
c wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >

c the integration is done along the imaginary axis
c w' ==> iw', w' is now real
c wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >

c transform: x = 1/(1+w')
c wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >

c the integrand is peak around w'=0 or x=1 when we=w-e=0
c to handel the problem, add and substract the singular part as follows:
c wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
c                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
c        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c the second term of the integral can be done analytically, which
c results in the last term

c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2

c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points

      implicit real*8 (a-h,o-z)
      real(8):: x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum,img=(0d0,1d0),wintz
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/

c if w = e the integral = -v(0)/2
cmod 93.06.23
c     if (dabs(we) < tol) then
c     wintz       = -0.5d0*v0
c     return
c     endif
      if (dabs(we) < tol) then
        print *, ' we=',we
Cstop2rx 2013.08.09 kino        stop 'wintz: |w-e| < tol'
        call rx( 'wintz: |w-e| < tol')
      endif
c frequency integral
      we2        = we*we
      sum        = 0.d0

      if (dabs(we) < rmax/a) then
        do       i = 1,nx
          omd   = 1d0/x(i) - 1d0
          onemx      = 1.d0 - x(i)
          cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + we*cons*(v(i) - v0*expa(i))*wt(i)
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)
          endif
        enddo

        wintz      = -sum/pi
     .            - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)
     .                                     *derfc(a*dabs(we))
      else
        do       i = 1,nx
          omd   = 1d0/x(i) - 1d0 !this was missing. I added this at 25June2008
          onemx      = 1.d0 - x(i)
          cons       = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + we*cons*v(i)*wt(i)
          if(npm==2) then !Asymmetric contribution for
            sum  = sum - cons* v(i+nx)*omd*wt(i)
          endif
        enddo

        wintz      = -sum/pi
      endif

      wintz_npm= wintz
      return
      end

c-------------------------------------------------------------
      complex(8) function wintz(v,v0,x,wt,a,expa,we,
     d                                nx)
c takao complex version of wint by ferdi
c 92.02.24
c wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >

c the integration is done along the imaginary axis
c w' ==> iw', w' is now real
c wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >

c transform: x = 1/(1+w')
c wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >

c the integrand is peak around w'=0 or x=1 when we=w-e=0
c to handel the problem, add and substract the singular part as follows:
c wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
c                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
c        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c the second term of the integral can be done analytically, which
c results in the last term

c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2

c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points

      implicit real*8 (a-h,o-z)
      real(8):: x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/

c if w = e the integral = -v(0)/2
cmod 93.06.23
c     if (dabs(we) < tol) then
c     wintz       = -0.5d0*v0
c     return
c     endif
      if (dabs(we) < tol) then
        print *, ' we=',we
Cstop2rx 2013.08.09 kino        stop 'wintz: |w-e| < tol'
        call rx( 'wintz: |w-e| < tol')
      endif
c frequency integral
      we2        = we*we
      sum        = 0.d0

      if (dabs(we) < rmax/a) then
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*(v(i) - v0*expa(i))*wt(i)
        enddo

        wintz      = -sum/pi
     .           - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)
     .                                    *derfc(a*dabs(we))
      else
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*v(i)*wt(i)
        enddo

        wintz      = -sum/pi
      endif

      return
      end

c--------------------------------------------------------------------
      complex(8) function alagr3z_even (x,xi,fi)
      implicit none
c even function version of alagr3z
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      real(8)::  x,xi(3)
      complex(8) :: fi(3) ,alagr3z
      real(8) :: amatinv(3,3),amat(3,3),ratio,detxx
      integer i

c--- Assure that dimag(fi) is positive.
Cstop2rx 2013.08.09 kino      if(x<0d0) stop ' alagr3z: x<0d0'
      if(x<0d0) call rx( ' alagr3z: x<0d0')
      do i=1,3
        if(dimag(fi(i))>1d-8) then
          ratio= abs(dimag(fi(i))/dreal(fi(i)))
          if( ratio >1d-6) then !jan2004. I added a bit of allowance
            print *,' alagr3z_even: error! dimag(fi(i)) is positive',i,xi,fi
Cstop2rx 2013.08.09 kino            stop    ' alagr3z: error! dimag(fi(i)) is positive'
            call rx( ' alagr3z: error! dimag(fi(i)) is positive')
          endif
        endif
      enddo
c
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
C     call dinv33(amat,0,amatinv,detxx)
      call minv33(amat,amatinv)
      alagr3z = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )

c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
c
      if(dimag(alagr3z)>0d0 ) then
        alagr3z = dcmplx( dreal(alagr3z),0d0)
      endif
      alagr3z_even=alagr3z
      end

c
c-------------------------------------
      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,ntp0,nstate),zwz(ntp0,nstate)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,ntp0,nstate) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
      do  it = 1,nstate
        do itp = 1,ntp0
          zwz(itp,it) = sum( dconjg(zmel(1:ngb,itp,it))*CC(1:ngb,itp,it))
        enddo
      enddo
      deallocate(CC)
      end
c-------------------------------
      subroutine gen_ua(ratio,niw,freqx, expa_,ua_)
      implicit none
      integer(4):: ix,niw,verbose
      real(8):: ratio,freqw, freqx(niw), expa_(niw),ua_  !, ua_min=1d0
      ix=niw
      freqw = (1d0 - freqx(ix))/ freqx(ix)

c      if( log(ratio) > -1d0/freqw/ua_min**2) then
c         print *,' ratio =',ratio
c         print 'gen_ua: ratio >1 strange ratio...'
c         stop  'gen_ua: ratio >1 strange ratio...'
c        ua_ = ua_min
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_ua: ratio >1 strange ratio...'
Cstop2rx 2013.08.09 kino        stop  ' gen_ua: ratio >1 strange ratio...'
        call rx( ' gen_ua: ratio >1 strange ratio...')
      else
        ua_ = sqrt(-1d0/freqw*log(ratio))
      endif
      do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
      enddo
      end
      subroutine gen_uaa(ratio,freqx, ua_)
      implicit none
      integer(4):: ix,niw,verbose
      real(8):: ratio,freqw, freqx, ua_
      freqw = (1d0 - freqx)/ freqx
      if( ratio > 1d0 ) then
        print *, ' ratio =',ratio
        print *, ' gen_uaa: ratio >1 strange ratio...'
Cstop2rx 2013.08.09 kino        stop  ' gen_uaa: ratio >1 strange ratio...'
        call rx( ' gen_uaa: ratio >1 strange ratio...')
      else
        ua_ = sqrt(-1d0/freqw*log(ratio))
      endif
      end

      subroutine gen_expa(niw,freqx,ua_,  expa_)
      implicit none
      integer(4):: ix,niw
      real(8)::  freqw, freqx(niw), expa_(niw), ua_
      do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
      enddo
      end

c-------------------------------------------------------------------
      subroutine matzwz2(iSigma_en ,zw,zmel, ntq,nstate,ngb, zwz) !sf 23May02
C     use m_mpi, only: mpi__rank
      implicit none
      integer(4) :: nstate,ntq,itp,itpp,it,ngb,iSigma_en
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntq)
      complex(8) :: zwz(nstate,ntq,ntq) ,zdotc                   !sf 23May02
      complex(8), allocatable :: CC(:,:,:)
      complex(8),allocatable :: z1r(:,:),z2r(:,:),zwzi(:,:)
      integer(4)::verbose
      integer(4):: ivc
      external verbose

      allocate(CC(ngb,nstate,ntq) )
      !write(*,'("10.1 matzwz begin ",$)')
      !call cputid(0)
C     write(6,'(i3,a,10i7)')mpi__rank,'matzwz2 enter',ngb,ngb,ngb,nstate*ntq,ngb,nstate*ntq
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
      if(verbose()>=50) call cputid2(' complete W z product',0)
      !write(*,'("10.2 matzwz after matm ",$)')
      !call cputid(0)
      if (iSigma_en==1.or.iSigma_en==5) then                                       !sf 23May02
        zwz=0d0
        do itp = 1,ntq
          do  it = 1,nstate
c        zwz(it,itp,itp)=sum(dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
            zwz(it,itp,itp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itp),1 )
          enddo
        enddo
      elseif (iSigma_en==2 .or. iSigma_en==3 ) then
C       Original inline zwz_t,t',t'' = sum_c z_c,t,t' CC_c,t,t''
C        do itpp= 1,ntq
C          do itp = 1,ntq
C            do  it = 1,nstate
C              zwz(it,itp,itpp)=
C     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1 )
Cc     &    sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itpp) )      !sf 23May02
C            enddo
C          enddo
C        enddo
C        if(verbose()>=50) call cputid2(' complete z+ (W z) product',0)
CC       call zprm('inline zwz',zwz,nstate,nstate,ntq**2)

C       Vectorized version of zwz_t,t',t'' = sum_c zmel*_c,t,t' CC_c,t,t'
C       Write as standard matrix product of t-dependent arrays
C       zwzi(t)_t',t'' = zwz_t,t',t''
C       z1(t)_t',c     = zmel*_c,t,t'
C       z2(t)_c,t''    = CC_c,t,t'
C       zwzi(t)_t',t'' = sum_c z1(t)_t',c z2(t)_c,t''
C       zwzi(t)        = z1(t) z2(t)
        allocate(z1r(ntq,ngb),z2r(ngb,ntq),zwzi(ntq,ntq))
        do  it = 1, nstate
          do  itp = 1, ntq
            do  ivc = 1, ngb
              z1r(itp,ivc) = dconjg(zmel(ivc,it,itp))
              z2r(ivc,itp) = CC(ivc,it,itp)
            enddo
          enddo
C         call zprm('z1r',z1r,ntq,ntq,ngb)
C         call zprm('z2r',z2r,ngb,ngb,ntq)
          call zgemm('N','N',ntq,ntq,ngb,(1d0,0d0),z1r,ntq,
     .      z2r,ngb,(0d0,0d0),zwzi,ntq)
          do  itp = 1, ntq
            do itpp = 1, ntq
              zwz(it,itp,itpp) = zwzi(itp,itpp)
            enddo
          enddo
        enddo
        deallocate(z1r,z2r,zwzi)
        if(verbose()>=50) call cputid2(' complete z+ (W z) product',0)
C       call zprm('vectorized zwz',zwz,nstate,nstate,ntq**2)

      else
Cstop2rx 2013.08.09 kino        stop "sxcf, matzwz2: iSigma_en /= 0,1,2,3"
        call rx( "sxcf, matzwz2: iSigma_en /= 0,1,2,3")
      endif
      !write(*,'("10.3 matzwz end ",$)')

      deallocate(CC)
      end

      subroutine matzwz3(iSigma_en,npm,zw,zmel,ndz,ntq,nstate,ngb, zwz)
C- Matrix element <zmel | W | zmel>
C ----------------------------------------------------------------------
Ci Inputs
Ci   iSigma_en
Ci         : mode determining form of static approximation to SE
Ci         : 0 SE_nn'(ef)+imag integral:delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ci         : 1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ci         : 2 SE_nn'((e_n+e_n')/2)
Ci         : 3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
Ci         : 4 SE_nn'(ef)
Ci         : 5 delta_nn' SE_nn(e_n)
Ci   zw    :Screened Coulomb interaction W
Ci   zmel  :matrix elements of eigenfunctions and the product basis
Ci   ndz   :leading dimension of zw
Ci   ntq   :number of states n or n' where SE_nn' is to be made
Ci   nstate:number of intermediate states t in zmel, to be summed over
Ci   ngb   :number of product basis functions
Co Outputs
Co   zwz   :<zmel | W | zmel>
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   matzwz2 was rewritten to make optimimum use of the BLAS3 library.
Cr   Can it be rewritten to avoid the allocation of CC ?
Cu Updates
Cu   1 May 13 Adapted from matzwz2
C ----------------------------------------------------------------------
      use fcugw
      use omp_lib
#ifdef MKL
      use mkl_service
#endif
      implicit none
C ... Passed parameters
      integer(4) :: ndz,nstate,ntq,ngb,iSigma_en,npm
      complex(8) :: zw(ndz,ngb),zmel(ngb,ntq,nstate)
      complex(8) :: zwz(ntq,ntq,nstate)

C ... Local parameters
      integer itp,itpp,it,ig, mkl_num_threads, ompid, ompnt(2), omp_gemm
      integer(4)::verbose,iverbose,idalloc
!       complex(8), allocatable :: CC(:,:,:)
      complex(8), allocatable :: cc(:,:), cc1(:,:,:)
      complex(8), allocatable :: z1r(:,:),z2r(:,:),zwzi(:,:)
      complex(8) :: zdotc
C ... External calls
      external cputid2,matm,zaxpy,zgemm

      integer(8) :: cr,cm,c1,c2
      call system_clock(count_rate = cr, count_max = cm)
!       call system_clock(c1)
!       call fcuzwz(ngb, ntq, nstate, zw, ndz, zmel, ngb, zwz, ntq, 0)
!       call system_clock(c2)
!       write(*,'("zwz2:",x,g20.12,"s")') real(c2 - c1,8)/real(cr,8)
!       return

C     iverbose = (verbose()/50)*10
      iverbose = 0
      ig = idalloc('zmw',iverbose+2,ngb,ntq)
!       allocate(cc(ngb,ntq,nstate))
!      allocate(cc(ngb,ntq))

!       print *, 'isigma_en, ngb, nstate, ntq', isigma_en, ngb, nstate, ntq

!       call system_clock(c1)

C ... Make CC(I,t,t') = S[J] W(I,J) ME(J,t,t')
! assuming zw is hermitian, then 'c' is conceptually the same as 'n' but is faster
!       call zgemm('c','N',ngb,nstate*ntq,ngb,(1d0,0d0),zw,ndz,zmel,ngb,(0d0,0d0),cc,ngb)
!       call system_clock(c2)

!       write(*,'("zwz1:",x,g20.12,"s")') real(c2 - c1,8)/real(cr,8)


      if(verbose()>=50) call cputid2(' complete W z product',0)

!       call system_clock(c1)

      if (iSigma_en==1 .or. iSigma_en==5) then !sf 23May02
        allocate(cc(ngb,ntq))
        zwz = 0d0
        do  it = 1, nstate
          call zgemm('c','n',ngb,ntq,ngb,(1d0,0d0),zw,ndz,zmel(1,1,it),ngb,(0d0,0d0),cc,ngb)
          do  itp = 1, ntq
c           zwz(it,itp,itp)=sum(dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
            zwz(itp,itp,it) = zdotc(ngb,zmel(1,itp,it),1,cc(1,itp),1)
          enddo
        enddo
        deallocate(cc)
      elseif (iSigma_en==2 .or. iSigma_en==3 ) then
C       Original inline zwz_t,t',t'' = sum_g z_g,t,t' CC_g,t,t''
C        do  itpp = 1, ntq
C          do  itp = 1, ntq
C            do  it = 1, nstate
C              zwz(it,itp,itpp)=
C     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1 )
Cc     &    sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itpp) )      !sf 23May02
C            enddo
C          enddo
C        enddo
C        if(verbose()>=50) call cputid2(' complete z+ (W z) product',0)
CC       call zprm('inline zwz',zwz,nstate,nstate,ntq**2)

C       Vectorized version of zwz_t,t',t'' = S[I] zmel*_I,t,t' CC_I,t,t''
C       Write as standard matrix product of t-dependent arrays
C       zwzi(t)_t',t'' = zwz_t,t',t''
C       z1(t)_t',g     = zmel*_g,t,t'
C       z2(t)_g,t'     = CC_g,t,t'
C       zwzi(t)_t',t'' = sum_g z1(t)_t',g z2(t)_g,t''
C       In matrix form:
C       zwzi(t)        = z1(t) z2(t)
!        allocate(z1r(ntq,ngb),z2r(ngb,ntq),zwzi(ntq,ntq))
!         allocate(zwzi(ntq,ntq))

        omp_gemm = 1
! Would have been nice to be able to use OMP_NUM_THREADS=N,4 but other places need the full N*4
#ifdef MKL
        omp_gemm = 4
#endif
        ompnt = [omp_get_max_threads(),omp_gemm]

        ompnt(2) = min(ompnt(1),ompnt(2))
        if (mod(ompnt(1),ompnt(2)) /= 0) stop 'number of threads not multiple of thread subteam (4).'
        ompnt(1) = ompnt(1)/ompnt(2)

        allocate(cc1(ngb,ntq,0:ompnt(1)-1))
        call omp_set_nested(.true.)

#ifdef MKL
        mkl_num_threads = mkl_get_max_threads()
        call mkl_set_dynamic(0)
        call mkl_set_num_threads(ompnt(2))
#endif

!$omp parallel num_threads(ompnt(1)) private(ompid)
        ompid = omp_get_thread_num()
!$omp do
        do  it = 1, nstate
          if (npm == 1) then
            call zhemm('l','l',ngb,ntq,(1d0,0d0),zw,ndz,zmel(1,1,it),ngb,(0d0,0d0),cc1(1,1,ompid),ngb)
          else
! the straightforward implementation would use zgemm type 'n,n' but the 'c,n' used here is compensated later in mkseci
            call zgemm('c','n',ngb,ntq,ngb,(1d0,0d0),zw,ndz,zmel(1,1,it),ngb,(0d0,0d0),cc1(1,1,ompid),ngb)
          endif
#ifdef Z3M
          call zgemm3m('c','n',ntq,ntq,ngb,(1d0,0d0),zmel(1,1,it),ngb,cc1(1,1,ompid),ngb,(0d0,0d0),zwz(1,1,it),ntq)
#else
          call zgemm('c','n',ntq,ntq,ngb,(1d0,0d0),zmel(1,1,it),ngb,cc1(1,1,ompid),ngb,(0d0,0d0),zwz(1,1,it),ntq)
#endif
        enddo
!$omp end do
!$omp end parallel
        deallocate(cc1)
        call omp_set_nested(.false.)

#ifdef MKL
        call mkl_set_dynamic(1)
        call mkl_set_num_threads(mkl_num_threads)
#endif

!        deallocate(z1r,z2r,zwzi)
!         deallocate(zwzi)
        if(verbose()>=50) call cputid2(' complete z+ (W z) product',0)
C       call zprm('vectorized zwz',zwz,nstate,nstate,ntq**2)

      else
        stop "sxcf, matzwz3: iSigma_en /= 1,2,3,5"
      endif

!       call system_clock(c2)
!       write(*,'("zwz2:",x,g20.12,"s")') real(c2 - c1,8)/real(cr,8)

!      deallocate(cc)
      ig = idalloc('zmw',iverbose+4,ngb,ntq)

      end

c--------------------------------------------------------------------
      complex(8) function alagr3z2(x,xi,fi,i,j)  !sf 24 May02
      implicit none
C     implicit real*8 (a-h,o-z)
c even function version of alagr3z
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      integer(4) :: i,j,n
      real(8)::  x,xi(3)
      complex(8) :: fi(3)
      real(8) :: amatinv(3,3),amat(3,3),ratio,detxx

c--- Assure that dimag(fi) is positive.
Cstop2rx 2013.08.09 kino      if(x<0d0) stop ' alagr3z2: x<0d0'
      if(x<0d0) call rx( ' alagr3z2: x<0d0')
      if (i==j) then  !sf 24 May02 only diagonal Im[Sigma_ii] should be negative
        do n=1,3
          if(dimag(fi(n))>1d-8) then
C           ratio= abs(dimag(fi(i))/dreal(fi(i)))
            ratio= abs(dimag(fi(n))/dreal(fi(n))) ! Bug fix MvS May 2013
            if( ratio >1d-8) then !jan2004. I added a bit of allowance
              print *,' alagr3z2: error! dimag(fi(i)) is positive',i,xi,fi
Cstop2rx 2013.08.09 kino              stop    ' alagr3z: error! dimag(fi(i)) is positive'
              call rx( ' alagr3z: error! dimag(fi(i)) is positive')
            endif
          endif
        enddo
      endif !(i=j)
c
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
C     call dinv33(amat,0,amatinv,detxx)
      call minv33(amat,amatinv)
      alagr3z2 = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )

c--- imag part of alagr3z may be positive because interpolation could yield
c    positive value even if all fi(i) <0.
c
      if(i==j) then
        if(dimag(alagr3z2)>0d0 ) then
          alagr3z2 = dcmplx( dreal(alagr3z2),0d0)
        endif
      endif !(i==j)
      end

      subroutine timeshow(info)
      character*(*) :: info
      write(*,'(a,$)')info
      call cputid(0)
      end




c======================================================================
      subroutine matzwz_notusedx(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
      do itp = 1,ntp0
        do  it = 1,nstate
          zwz(it,itp) = sum(
     &   dconjg(zmel(:,it,itp))*matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)))
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
        enddo
      enddo
      end
c---------------------------------------------
      subroutine matzwzs_notused(zw,zmel, ntq,nstate,ngb, zwz) !sf 23May02
c I have not checked it yet...
      implicit none
      integer(4) :: nstate,ntq,itp,itpp,it,ngb,iSigma_en,igb,igb2
      complex(8) :: zw(ngb,ngb),zwx(ngb,ngb),zmel(ngb,nstate,ntq)
      complex(8) :: zwz(nstate,ntq,ntq) ,zdotc,zz3(ntq),zz2
      complex(8), allocatable :: CC(:,:,:),zwzx(:,:,:),zmel1_(:,:,:)
      allocate(zwzx(nstate,ntq,ntq),zmel1_(ntq,ngb,nstate))
      do itpp= 1,ntq
        do it  = 1,nstate
          zmel1_(itpp,1:ngb,it) = zmel(1:ngb,it,itpp)
        enddo
      enddo
      zwx = zw
      do igb=1,ngb
        zwx(igb,igb)=.5d0*zwx(igb,igb)
      enddo
      do itp=1,ntq
        do it =1,nstate
c        zmel1(:)=dconjg(zmel(:,it,itp))
          zz3=0d0
          do igb2=1,ngb
            zz2 = zdotc(igb2, zmel(1,it,itp),1,zwx(1,igb2),1 )
            call zaxpy(ntq,zz2,zmel1_(1,igb2,it),1,zz3,1)
          enddo
          zwzx(it,itp,:)= zz3(:)
        enddo
      enddo
      do itpp=1,ntq
        do itp =1,ntq
          do it  =1,nstate
            zwz(it,itp,itpp)=zwzx(it,itp,itpp)+dconjg(zwzx(it,itpp,itp))
          enddo
        enddo
      enddo
      deallocate(zwzx,zmel1_)
      return
c-----------------------------------------------
      Allocate(CC(ngb,nstate,ntq) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq) !most time consuming part
      !write(*,'("10.2 matzwz after matm ",$)')
      !call cputid(0)
      do itpp= 1,ntq
        do itp = 1,ntq
          do  it = 1,nstate
            zwz(it,itp,itpp)=
     &    zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1 )
          enddo
        enddo
      enddo
      deallocate(CC)
      end
c      subroutine matzwz(zw,zmel, ntq,nstate,ngb, zwz)
c      implicit none
c      integer(4) :: nstate,ntq,itp,it,ngb
c      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntq),zwz(nstate,ntq)
c      complex(8), allocatable :: CC(:,:,:)
c      allocate(CC(ngb,nstate,ntq) )
c
c      call matm(zw,zmel,cc, ngb, ngb, nstate*ntq)
c
c      do itp = 1,ntq
c      do  it = 1,nstate
c        zwz(it,itp) = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
c      enddo
c      enddo
c
c      deallocate(CC)
c      end

c-------------------------------------
c   This does not work.
c      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
c      implicit none
c      integer(4) :: nstate,ntp0,itp,it,ngb
c      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
c        do itp = 1,ntp0
c        do  it = 1,nstate
c        call mml(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp), wx(1:ngb))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp))*wx )
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
c        enddo
c        enddo
c      end

c -----------------------------------
      complex(8) function alagr3z (x,xi,fi)
c complex version of fi
c 92.03.02
c three-point interpolation with arbitrary mesh
c given a function fi at three points xi, the routine interpolates
c the function at x
c f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
c      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
c      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3
c
c x  = the point at which the function is to be interpolated
c xi(3) = points where the function is given
c fi(3) = the function at xi
      implicit real*8 (a-h,o-z)
      real(8):: xi(3)
      complex(8) :: fi(3)
      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x23        = xi(2)-xi(3)
      alagr3z    = xx2*xx3*fi(1)/(x12*x13)
     .           - xx1*xx3*fi(2)/(x12*x23)
     .           + xx1*xx2*fi(3)/(x13*x23)
      return
      end

c--------------------------------------------------------------------
      complex(8) function alagr3zz(x,xi,fi)  !jan 2008 takao. starting from alagr3z2
c evne function version of alagr3z !
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      implicit none
      real(8)::  xi(3), amatinv(3,3),amat(3,3),x,detxx
      complex(8) :: fi(3)
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
C     call dinv33(amat,0,amatinv,detxx)
      call minv33(amat,amatinv)
      alagr3zz = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )
c      write(6,"(' coeff amatinv = ',9d11.3)") matmul(amatinv,dimag(fi))
c      write(6,"(' coeff    x    = ',9d11.3)") x**2,x**4
c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
      if(dimag(alagr3zz)>0d0) alagr3zz = dcmplx( dreal(alagr3zz),0d0)
      end

      complex(8) function wintz0(npm,ix,v,v0,x,wt,a,expa,we,nx)
C- Integral [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) done along imaginary axis
C ----------------------------------------------------------------------
Ci Inputs
Ci   npm   :1 Use time reversal symmetry
Ci         :2 No time reversal symmetry
Ci   ix    :0 Return entire integral
Ci         :>0 Return contribution to integral from point ix
Ci   v     :Function tabulated on nx points iw' on the imaginary axis; see Remarks
Ci   v0    :v(omega=0)
Ci   x     :Vector of transformed frequencies; see Remarks
Ci   wt    :Integration weights
Ci   a     :Inverse width of gaussian to be integrated out analytically,
Ci         :making residual small; see Remarks
Ci   expa  :Vector of exp(-a^2 w'^2) for the nx frequencies iw'
Ci   we    :w - e
Ci   nx    :Number of points in quadrature
Co Outputs
Co  wintz0 :Integral
Cl Local variables
Cl         :
Cr Remarks
Cr   The integration is done along the imaginary axis
Cr   w' ==> iw', w' is now real.
Cr   wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >
Cr
Cr   Transform: x = 1/(1+w')
Cr     wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
Cr   This results in a denser mesh in w' around 0 for quadrature in x
Cr   Useful when the integrand is peaked around w'=0 or x=1 when we=w-e=0.
Cr   To make the integration smooth, add and substract the singular part as follows:
Cr     wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
Cr                        *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
Cr             - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
Cr   The second term of the integral can be done analytically, which
Cr   results in the last term.
Cr
Cr   When w-e -> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} -> delta(w')
Cr   The integral reduces to -v(0)/2
Cu Updates
Cu   10 May 13 Adapted from Takao's wintz_npm
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ix,npm,nx
      real(8):: we,a,x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx,2),v0
C ... Local parameters
      integer i,low,high
      real(8) pi,tol,we2,omd,onemx,cons
      real(8) :: rmax=2d0 !rmax=2 chosen by Ferdi. Is it optimum? See wintz
      real(8) derfc
      external derfc
      complex(8) :: sum,wintz

      data pi/3.1415926535897932d0/, tol/1d-8/

      if (dabs(we) < tol) then
        wintz0 = -v0/2
C       Until we can check it
        print *,' we=',we
        stop 'wintz0: |w-e| < tol'
      endif

C ... Setup
      we2 = we*we
      sum = 0d0
      if (ix == 0) then
        low = 1; high = nx
      else
        low = ix; high = ix
      endif

      if (dabs(we) < rmax/a) then
        do  i = low, high
          omd = 1d0/x(i) - 1d0
          onemx = 1 - x(i)
          cons = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum = sum + we*cons*(v(i,1) - v0*expa(i))*wt(i)
C         Asymmetric contribution when no time reversal
          if (npm == 2) then
            sum = sum - cons*v(i,2)*omd*wt(i)
          endif
        enddo
        if (ix <= 1) then     ! ix-independent
          wintz = -v0/2*dsign(1d0,we)*dexp(we2*a*a)*derfc(a*dabs(we))
        else
          wintz = 0
        endif
        wintz = wintz - sum/pi

      else
        do  i = low, high
          omd = 1d0/x(i) - 1d0
          onemx = 1 - x(i)
          cons = 1d0/(we2*x(i)*x(i) + onemx*onemx)
          sum = sum + we*cons*v(i,1)*wt(i)
C         Asymmetric contribution when no time reversal
          if (npm==2) then
            sum = sum - cons* v(i,2)*omd*wt(i)
          endif
        enddo
        wintz = -sum/pi
      endif

      wintz0 = wintz
      end

      complex(8) function wintzg(npm,ix,v,v0,x,wt,a,expa,we,nx,esmr)
C- Gaussian smeared-pole version of wintz0.
C ----------------------------------------------------------------------
Ci Inputs
Ci   npm   :1 Use time reversal symmetry
Ci         :2 No time reversal symmetry
Ci   ix    :0 Return entire integral
Ci         :>0 Return contribution to integral from point ix
Ci   v     :Function tabulated on nx points iw' on the imaginary axis; see Remarks
Ci   v0    :v(omega=0)
Ci   x     :Vector of transformed frequencies; see Remarks
Ci   wt    :Integration weights
Ci   a     :constant determining the range of exp(-a^2 w'^2)
Ci         :for analytic part of integral; see Remarks
Ci   expa  :exp(-a^2 w'^2)
Ci   we    :w - e
Ci   nx    :Number of points in quadrature
Ci   esmr  :smearing of eigenvalue (Ry)
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   Assume that each eigenvalue is broadened by width of esmr
Cu Updates
Cu   10 May 13
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ix,npm,nx
      real(8):: we,a,esmr,x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx,2),v0
C ... Local parameters
      integer i,low,high
      real(8) pi,we2,omd,cons,sig2,omd2w,derfcx,sig,aw,eee,aw2
      complex(8) :: sum,wintz0,sumgauss
C     complex(8) :: wwz
C     integer :: ie,nav = 2000
      data pi/3.1415926535897932d0/

C ... Brute-force integration (for checking the fast version)
C      do  ie = 1, nav
C        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
C        www = 1d0
C        if (ie==1 .or. ie==nav) www = 0.5d0
C        wwz = wwz + www*wintz0(npm,ix,v,v0,x,wt,a,expa,wex,nx)
C      enddo
C      wintzg = wwz/(nav-1)
C      return

C ... No broadening of the pole: just wintz0
      if (esmr == 0d0) then
        wintzg = wintz0(npm,ix,v,v0,x,wt,a,expa,we,nx)
        return
      endif

C ... Setup
      sig  = .5d0*esmr
      sig2 = 2*sig**2
      we2  = we**2
      sum  = 0

C --- Numerical integration of nonsingular part on the imaginary axis ---
C     Pole weight is given by:
C     1-exp(-R^2/sig2) = 1/N \int_0^R exp(- x^2/sig2) 2 pi r dr
      if (ix == 0) then
        low = 1; high = nx
      else
        low = ix; high = ix
      endif
      do  i = low, high
        omd = 1d0/x(i) - 1d0
        omd2w = omd**2 + we2

C       Gauss theorem ... but is this correct? Not three dimensional...
        if (omd2w/sig2 > 5d-3) then
          cons = (1d0 - exp(-omd2w/sig2))/omd2w
        else
          cons = (1d0/sig2 - omd2w/sig2**2/2d0
     .            + omd2w**2/sig2**3/6d0 - omd2w**3/sig2**4/24d0
     .            + omd2w**4/sig2**5/120d0 - omd2w**5/sig2**6/720d0)
        endif
        sum  = sum + we*cons * (v(i,1)-v0*expa(i)) * wt(i)/(x(i)**2)
C       Asymmetric contribution when no time reversal
        if (npm==2) then
          sum = sum - cons*v(i,2)*omd*wt(i)/(x(i)**2)
        endif
      enddo


C --- Singular part ---
      if (ix <= 1) then ! ix-independent
        aw = abs(a*we)
        aw2 = aw**2
        if (aw < 6d0) then
          eee = we**2/sig2
          sumgauss = dsign(1d0,we)*.5d0*v0*exp(aw2)*
     .      (derfcx(sqrt(aw2+eee)) - derfcx(aw) )
        else !April2004
          sumgauss = -.5d0*v0*(2d0/(a*we)/pi)*
     .      dsqrt(pi)/2d0*
     .       (1d0 - 1d0/2d0/aw2*
     .         (1d0 - 3d0/2d0/aw2*
     .           (1d0 - 5d0/2d0/aw2*
     .             (1d0 - 7d0/2d0/aw2*
     .               (1d0 - 9d0/2d0/aw2*
     .                 (1d0 - 11d0/2d0/aw2
     .     ))))))
        endif
      else
        sumgauss = 0
      endif
      wintzg = -sum/pi + sumgauss

C      test code
C      print *,' sum= ', sum
C      print *,' awe1= ', a*we
C      print *,' erf1= ', derfcx(a*we)
C      print *,' awe2= ', sqrt(a**2+1d0/sig2)*we
C      print *,' erf2= ', derfcx(sqrt(a**2+1d0/sig2)*we)
C      print *,' dexp= ', dexp(we2*a**2)
C      print *,' erf3= ', derfcx(aw)
C      print *,' erf4= ', derfcx(sqrt(aw2+eee))
C      print *,' aw2= ', aw2
C      aw = a*we
C      write(*, "(' wintzg chk1=',2f8.3, 4d18.10)" )we,sig,
C     . -.5d0*v0*dexp(aw**2)*( derfcx(aw)-derfcx(sqrt(a**2+1d0/sig2)*we))
C      write(*, "(' wintzg chk2=',2f8.3, 4d18.10)" )we,sig,sumgauss1
C      write(*, "(' wintzg chk3=',2f8.3, 4d18.10)" )we,sig,sumgauss2
      end
C      subroutine fmain
C      implicit none
C      integer ix,npm,nx
C      parameter (nx=6)
C      real(8):: we,a,esmr,x(nx),wt(nx),expa(nx),freqx(nx)
C      complex(8) ::v(nx,2),v0
C      real(8):: tol=1d-12
C      double precision freqw
C      complex(8) :: wres,wres2,wintzg
C
C      v(1,1) = (-9.8436671501287847d-05, 0d0)
C      v(2,1) = (-0.0019192082020032866d0, 0d0)
C      v(3,1) = (-0.0047310624388582296d0, 0d0)
C      v(4,1) = (-0.0076706731521482404d0, 0d0)
C      v(5,1) = (-0.014478269823812832d0, 0d0)
C      v(6,1) = (-0.026873663652901166d0, 0d0)
C      v0 = (-0.031138884443444498d0, 0d0)
C      we = -0.0014037973777867307d0
C      esmr = .003d0
C      npm = 1
C      x(1) = 0.033765242898423975d0
C      x(2) = 0.16939530676686748d0
C      x(3) = 0.38069040695840151d0
C      x(4) = 0.61930959304159849d0
C      x(5) = 0.83060469323313257d0
C      x(6) = 0.96623475710157603d0
C      wt(1) = 0.085662246189584998d0
C      wt(2) = 0.1803807865240695d0
C      wt(3) = 0.23395696728634549d0
C      wt(4) = 0.23395696728634549d0
C      wt(5) = 0.1803807865240695d0
C      wt(6) = 0.085662246189584998d0
C      freqx(1) = 0.033765242898423975d0
C      freqx(2) = 0.16939530676686748d0
C      freqx(3) = 0.38069040695840151d0
C      freqx(4) = 0.61930959304159849d0
C      freqx(5) = 0.83060469323313257d0
C      freqx(6) = 0.96623475710157603d0
C      a = 1
C      do  ix = 1, nx
C        freqw     = (1d0 - freqx(ix))/ freqx(ix)
C        expa(ix) = exp(-(a*freqw)**2)
C      enddo
C
C      wres = wintzg(npm,0,v,v0,x,wt,a,expa,we,nx,esmr)
C      call info0(0,0,-1,'test 1 ... compare against known result ...')
C      if (abs(wres-(-0.009911871046332333d0, 0d0)) < tol) then
C        print *, 'ok'
C      else
C        call rx('failed test')
C      endif
C
C      call info0(0,0,-1,
C     .  'test 2 ... accumulate integral point by point ...')
C      wres = 0
C      do  ix = 1, nx
C        wres = wres + wintzg(npm,ix,v,v0,x,wt,a,expa,we,nx,esmr)
C      enddo
C      if (abs(wres-(-0.009911871046332333d0, 0d0)) < tol) then
C        print *, 'ok'
C      else
C        call rx('failed test')
C      endif
C
C      call info0(0,0,-1,'test 3 ... check esmr=0 ...')
C      wres = 0
C      do  ix = 1, nx
C        wres = wres + wintzg(npm,ix,v,v0,x,wt,a,expa,we,nx,0d0)
C      enddo
C      if (abs(wres-wintzg(npm,0,v,v0,x,wt,a,expa,we,nx,0d0)) < tol)then
C        print *, 'ok'
C      else
C        call rx('failed test')
C      endif
C
C      call info0(0,0,-1,'test 4 ... esmr=0, below w=e>rmax/a ...')
C      a = 2d0-1e-6
C      do  ix = 1, nx
C        freqw     = (1d0 - freqx(ix))/ freqx(ix)
C        expa(ix) = exp(-(a*freqw)**2)
C      enddo
C      we = 1
C      wres = 0
C      do  ix = 1, nx
C        wres = wres + wintzg(npm,ix,v,v0,x,wt,a,expa,we,nx,0d0)
C      enddo
C      if (abs(wres-wintzg(npm,0,v,v0,x,wt,a,expa,we,nx,0d0)) < tol)then
C        print *, 'ok'
C      else
C        call rx('failed test')
C      endif
C      wres2 = wres
C
C      call info0(0,0,-1,'test 5 ... esmr=0, above w=e>rmax/a ...')
C      a = 2d0+1e-6
C      do  ix = 1, nx
C        freqw     = (1d0 - freqx(ix))/ freqx(ix)
C        expa(ix) = exp(-(a*freqw)**2)
C      enddo
C      we = 1
C      wres = 0
C      do  ix = 1, nx
C        wres = wres + wintzg(npm,ix,v,v0,x,wt,a,expa,we,nx,0d0)
C      enddo
C      if (abs(wres-wintzg(npm,0,v,v0,x,wt,a,expa,we,nx,0d0)) < tol)then
C        print *, 'ok'
C      else
C        call rx('failed test')
C      endif
C      print *, 'difference in the two:', wres2-wres
C
C
CC      do  esmr = 0, 1d-2, 1d-4
CC        wres = wintzg(npm,0,v,v0,x,wt,a,expa,we,nx,esmr)
CC        print *, esmr, dble(wres), dimag(wres)
CC      enddo
CC      stop
C
C      end
C
C      real(8) function derfcx(a)
C      real(8):: a, derfc,ax, amx=12d0
C      ax = a
C      if( abs(a)>amx) ax= a/abs(a) * amx
Cc      print *,' xxx ',ax
Cc      print *,' yyy ',derfc(ax)
C      derfcx=derfc(ax)
C      end
