      subroutine rdeq_old(enum1,enum2,ksop,z,v,rofi,nr,nsp,a,b,l,lmx,
     .                imu,scalede,gmt,fmt,gmtde,fmtde,gsmt,pprel)
C-Solves radial Dirac equations, in the presence of an effective magnetic field
C ----------------------------------------------------------------------
Ci Inputs
Ci   enum1 :linearization energy, spin 1
Ci   enum2 :linearization energy, spin 2
Ci   ksop  :spin orbit coupling parameters, used to adjust enu1,enu2
Ci   z     :nuclear charge
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l     :l quantum number
Ci   lmx   :dimensions ksop and pprel
Ci   imu   :projection of j (positive integer, counts -l-1/2...l+1/2)
Ci   scalede: numerical factor for energy differentiation
Co Outputs
Co   gmt   :(2x2) g amplitudes x r at rmax
Co   fmt   :(2x2) f amplitudes x r at rmax
Co   gsmt  : ?
Co   pprel :Relativistic potential parameters in kappa-mu representation
Co         :pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)
Co         :Not made here:
Co         :pprel(1,:,:,1:2,1:2) = C
Co         :pprel(2,:,:,1:2,1:2) = gamma
Co         :pprel(3,:,:,1:2,1:2) = D+, where (D+) D = delta
Co         :Made here:
Co         :pprel(4,:,:,1:2,1:2) = small parameter p
Co   gmtde :(2x2 matrix) gdot at rmax
Co   fmtde :(2x2 matrix) fdot at rmax
Cl Local variables
Cr Remarks
Cr   In this version enu depends only on kappa
Cu Updates
Cu   18 Oct 13 (Belashchenko) cleanup, bug fix
Cu   18 Jun 04 (A Chantis) working version
Cu   24 Jun 03 (A Chantis) adopted from V. Antropov
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,imu,l,lmx
      double precision v(nr,nsp),rofi(nr),z,a,b,enum1,enum2,scalede
      double precision ksop(0:lmx,nsp,nsp,6)
      double precision pprel(5,0:lmx,2*(lmx+1),2,2)
      double precision gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
C ... Local parameters
      integer n,np,ie,alpha,a1,nrk,it,nrmx
C     integer k1,k2,i,alp,ifi,fopn
C     character*(10) outs*80, fil*12
      double precision J,    !Jacobian
     .  c,csq,               !Speed of light 274.072d0
     .  de,                  !Finite energy difference for num diff.
     .  bmint,df11,df12,df13,df14,df21,df22,df23,df24,dg11,dg12,dg13,
     .  dg14,dg21,dg22,dg23,dg24,dx,e1,e2,eav,enu1,enu2,f1c,f1p,f1sn,
     .  f1snf1,f1snm1,f1snm2,f1snm3,f1snm4,f2c,f2p,f2sn,f2snf2,f2snm1,
     .  f2snm2,f2snm3,f2snm4,ff,g1c,g1p,g1sn,g1sng1,g1snm1,g1snm2,
     .  g1snm3,g1snm4,g2c,g2p,g2sn,g2sng2,g2snm1,g2snm2,g2snm3,g2snm4,
     .  gamma,hoc,mu,norm,nuf1int,nuf2int,nug1int,nug2int,pr,prodr,q,r,
     .  r1,rmax,smalpa(2,2),sqru,u,up,upp,vmint,w0,w1,w2,w3,x1,x2,xk,xk1
      real(8),parameter :: NULLR=-99999d0
      double precision g1(nr),g2(nr),f1(nr),f2(nr),g1s(nr),g2s(nr),
     .  f1s(nr),f2s(nr),nug1(nr),nug2(nr),nuf1(nr),nuf2(nr),bm(nr),
     .  vm(nr),psi(nr,4,2,5),psd(nr,4,2)
C ... External calls
      external dinv22,mul22,rx
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      eav = (enum1 + enum2)/2
      enu1 = eav
      enu2 = eav
      if (ksop(0,1,1,1) /= NULLR) then
        enu1 = -(l+1) * ksop(l,1,1,1)/2 + eav
        enu2 =   l    * ksop(l,1,1,1)/2 + eav
      endif

      mu  = imu - l - 1.5d0

C --- Test: Eliminate Spin-Orbit for non-mag. materials
C     mu = dble(l)+0.5d0

      if (z == 0) return

      csq = c*c
      nrk = 6 ; nrmx = nr
      hoc = (2*z/c)
      de = 0.03d0*scalede ; dx = 1d0
      ff = 1

      rmax = rofi(nr)
      vm(:) = (v(:,1) + v(:,2))/2
      bm(:) = ((v(:,2) - v(:,1))/2)

      if (imu == 1 .or. imu == 2*l+2) then
        xk = -l-1
      else
        xk = l
      endif
      xk1 = -l-1

      up = mu/(l-0.5d0) ; u = mu/(l+0.5d0) ; upp = mu/(l+1.5d0)
      sqru = dsqrt(1-u*u)

C ... Added to take care of the abs(mu)=l+0.5 case
      if (imu == 1 .or. imu == 2*l+2) up = upp

      g1(1) = 0 ; f1(1) = 0 ; g2(1) = 0 ; f2(1) = 0

c     de = 0.15
c     print *,'de=',de
      do  alpha = 1, 2

C --- Initial conditions at r-->0: V = V0 - 2Z/r, B = const
C     The wave function is a polynomial g0i = r^(gi-1)*Sum_over_v(g1v*r^v)

      if (alpha == 1) then
        gamma = dsqrt(xk*xk - hoc*hoc)
        g1(2) = 1d0 ; f1(2) = (xk + gamma)/hoc
        g2(2) = 0 ; f2(2) = 0
      else
        gamma = dsqrt(xk1*xk1 - hoc*hoc)
        g1(2) = 0 ; f1(2) = 0
        g2(2) = 1d0 ; f2(2) = (xk1 + gamma)/hoc
      endif

      do  ie = 1, 5
        e2 = enu2 + (ie - 3)*de ; e1 = enu1 + (ie - 3)*de

      if (imu == 1 .or. imu == 2*l+2) then
        nug1(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq
     .             + ff*(e2 - vm(2:nr) - upp*bm(2:nr))/csq
        nuf1(2:nr) = - (e2 + 2*z/rofi(2:nr) - vm(2:nr) - u  *bm(2:nr))
      else
        nug1(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq + ff*(e1 - vm(2:nr)
     .             + up*bm(2:nr))/csq
        nuf1(2:nr) = - (e1 + 2*z/rofi(2:nr) - vm(2:nr) + u *bm(2:nr))
      endif

      nug2(2:nr) = 1d0 + 2*z/rofi(2:nr)/csq
     .           + ff*(e2 - vm(2:nr) - upp*bm(2:nr))/csq
      nuf2(2:nr) = - (e2 + 2*z/rofi(2:nr) - vm(2:nr) - u  *bm(2:nr))

C     nug1(1) = nug1(2); nug2(1) = nug2(2)
C     nuf1(1) = nuf1(2); nuf2(1) = nuf2(2)
C     Extrapolate to the origin
      call ratint(rofi(2),nug1(2),4,rofi(1),nug1(1),w0)
      call ratint(rofi(2),nuf1(2),4,rofi(1),nuf1(1),w0)
      call ratint(rofi(2),nug2(2),4,rofi(1),nug2(1),w0)
      call ratint(rofi(2),nuf2(2),4,rofi(1),nuf2(1),w0)

C --- Runge-Kutta for the first NKR points ------------------------------
      J = a*(rofi(2) + b) ; q = dexp(a/2) ; r = rofi(2)

      g1s(2) = J*(c*nug1(2)*f1(2)                       - xk /r*g1(2))
      f1s(2) = J*( (nuf1(2)*g1(2) - sqru*bm(2)*g2(2))/c + xk /r*f1(2))
      g2s(2) = J*(c*nug2(2)*f2(2)                       - xk1/r*g2(2))
      f2s(2) = J*( (nuf2(2)*g2(2) - sqru*bm(2)*g1(2))/c + xk1/r*f2(2))

C     First point
      n = 2

C     Re-entry for subsequent points
    2 continue
      g1c = g1(n) ; f1c = f1(n) ; g2c = g2(n) ; f2c = f2(n)
      dg11 = dx*J* (c*nug1(n)*f1c                     - xk /r*g1c)
      df11 = dx*J* ( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk /r*f1c)
      dg21 = dx*J* (c*nug2(n)*f2c                     - xk1/r*g2c)
      df21 = dx*J* ( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

      g1c = g1(n) + dg11/2 ; f1c = f1(n) + df11/2
      g2c = g2(n) + dg21/2 ; f2c = f2(n) + df21/2

      J = J*q ; r = J/a - b

      vmint = (5*vm(n) + 2*vm(n+1) + vm(n+2))/8
      bmint = (5*bm(n) + 2*bm(n+1) + bm(n+2))/8

c     vmint = (6*vm(n) + 5*vm(n+1) - 3*vm(n-1))/8
c     bmint = (6*bm(n) + 5*bm(n+1) - 3*vm(n-1))/8

      if (imu == 1 .or. imu == 2*l+2) then
        nug1int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
        nuf1int =     - 2*z/r     - (e2 - vmint - u   * bmint)
      else
        nug1int = 1d0 + 2*z/r/csq + ff*(e1 - vmint + up * bmint)/csq
        nuf1int =     - 2*z/r     - (e1 - vmint + u  * bmint)
      endif

      nug2int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
      nuf2int =     - 2*z/r     - (e2 - vmint - u   * bmint)


      dg12 = dx*J*(c*nug1int*f1c                     - xk /r*g1c)
      df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + xk /r*f1c)
      dg22 = dx*J*(c*nug2int*f2c                     - xk1/r*g2c)
      df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

      g1c = g1(n) + dg12/2 ; f1c = f1(n) + df12/2
      g2c = g2(n) + dg22/2 ; f2c = f2(n) + df22/2

      dg13 = dx*J*(c*nug1int*f1c                     - xk /r*g1c)
      df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + xk /r*f1c)
      dg23 = dx*J*(c*nug2int*f2c                     - xk1/r*g2c)
      df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + xk1/r*f2c)

      g1c = g1(n) + dg13 ; f1c = f1(n) + df13
      g2c = g2(n) + dg23 ; f2c = f2(n) + df23

      n = n + 1

      J = J*q ; r = J/a - b

      dg14 = dx*J*(c*nug1(n)*f1c                     - xk /r*g1c)
      df14 = dx*J*( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + xk /r*f1c)
      dg24 = dx*J*(c*nug2(n)*f2c                     - xk1/r*g2c)
      df24 = dx*J*( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + xk1/r*f2c)

      g1(n) = g1(n-1) + (dg11 + 2*(dg12+dg13) + dg14)/6
      f1(n) = f1(n-1) + (df11 + 2*(df12+df13) + df14)/6
      g2(n) = g2(n-1) + (dg21 + 2*(dg22+dg23) + dg24)/6
      f2(n) = f2(n-1) + (df21 + 2*(df22+df23) + df24)/6

C ... Derivatives dP/dx and dQ/dx----------------------------

      g1s(n) = J*(c*nug1(n)*f1(n)                       - xk /r*g1(n))
      f1s(n) = J*( (nuf1(n)*g1(n) - sqru*bm(n)*g2(n))/c + xk /r*f1(n))
      g2s(n) = J*(c*nug2(n)*f2(n)                       - xk1/r*g2(n))
      f2s(n) = J*( (nuf2(n)*g2(n) - sqru*bm(n)*g1(n))/c + xk1/r*f2(n))

C --- The rest of the integration by Milne's method

      if (n < nrk) goto 2
      if (n == nrmx) goto 12

      g1sn = g1s(n) ; f1sn = f1s(n) ; g2sn = g2s(n) ; f2sn = f2s(n)

      g1snm1 = g1s(n-1) ; f1snm1 = f1s(n-1)
      g2snm1 = g2s(n-1) ; f2snm1 = f2s(n-1)

      g1snm2 = g1s(n-2) ; f1snm2 = f1s(n-2)
      g2snm2 = g2s(n-2) ; f2snm2 = f2s(n-2)

      g1snm3 = g1s(n-3) ; f1snm3 = f1s(n-3)
      g2snm3 = g2s(n-3) ; f2snm3 = f2s(n-3)

      g1snm4 = g1s(n-4) ; f1snm4 = f1s(n-4)
      g2snm4 = g2s(n-4) ; f2snm4 = f2s(n-4)

    4 J = J*q*q ; r = J/a - b

      g1p = g1(n-5)
     .    + (3*dx/10)*(11*g1sn-14*g1snm1+26*g1snm2-14*g1snm3+11*g1snm4)
      f1p = f1(n-5)
     .    + (3*dx/10)*(11*f1sn-14*f1snm1+26*f1snm2-14*f1snm3+11*f1snm4)
      g2p = g2(n-5)
     .    + (3*dx/10)*(11*g2sn-14*g2snm1+26*g2snm2-14*g2snm3+11*g2snm4)
      f2p = f2(n-5)
     .    + (3*dx/10)*(11*f2sn-14*f2snm1+26*f2snm2-14*f2snm3+11*f2snm4)

      do  it = 1, 100
        g1sng1 = J*(c*nug1(n+1)*f1p                       - xk/r*g1p)
        f1snf1 = J*( (nuf1(n+1)*g1p - sqru*bm(n+1)*g2p)/c + xk/r*f1p)
        g2sng2 = J*(c*nug2(n+1)*f2p                       - xk1/r*g2p)
        f2snf2 = J*( (nuf2(n+1)*g2p - sqru*bm(n+1)*g1p)/c + xk1/r*f2p)

        g1c = g1(n-3)
     .      + (2*dx/45)*(7*g1sng1+32*g1sn+12*g1snm1+32*g1snm2+7*g1snm3)
        f1c = f1(n-3)
     .      + (2*dx/45)*(7*f1snf1+32*f1sn+12*f1snm1+32*f1snm2+7*f1snm3)
        g2c = g2(n-3)
     .      + (2*dx/45)*(7*g2sng2+32*g2sn+12*g2snm1+32*g2snm2+7*g2snm3)
        f2c = f2(n-3)
     .      + (2*dx/45)*(7*f2snf2+32*f2sn+12*f2snm1+32*f2snm2+7*f2snm3)

        if (dabs(g1c-g1p) <= dabs(g1c)*1d-12 .and.
     .      dabs(f1c-f1p) <= dabs(f1c)*1d-12 .and.
     .      dabs(g2c-g2p) <= dabs(g2c)*1d-12 .and.
     .      dabs(f2c-f2p) <= dabs(f2c)*1d-12) exit

C ...   Check for convergence
        if (it == 100) then
          print *,abs(g1c-g1p),abs(f1c-f1p),abs(g2c-g2p),abs(f2c-f2p)
          call rx('Convergence not reached after 100 iterations')
        endif
        g1p = g1c ; f1p = f1c ; g2p = g2c ; f2p = f2c
      enddo

      n = n + 1

      g1(n) = g1c ; f1(n) = f1c ; g2(n) = g2c ; f2(n) = f2c
      g1s(n) = g1sng1; f1s(n) = f1snf1; g2s(n) = g2sng2; f2s(n) = f2snf2
      g1snm4 = g1snm3; f1snm4 = f1snm3; g2snm4 = g2snm3; f2snm4 = f2snm3
      g1snm3 = g1snm2; f1snm3 = f1snm2; g2snm3 = g2snm2; f2snm3 = f2snm2
      g1snm2 = g1snm1; f1snm2 = f1snm1; g2snm2 = g2snm1; f2snm2 = f2snm1
      g1snm1 = g1sn ; f1snm1 = f1sn ; g2snm1 = g2sn ; f2snm1 = f2sn
      g1sn = g1sng1 ; f1sn = f1snf1 ; g2sn = g2sng2 ; f2sn = f2snf2

      if (n < nrmx) goto 4

   12 continue

C ... Save and normalize the Psi functions
      psi(:,1,alpha,ie) = g1(:) ; psi(:,2,alpha,ie) = f1(:)
      psi(:,3,alpha,ie) = g2(:) ; psi(:,4,alpha,ie) = f2(:)
      pr = prodr(0,psi(1,1,alpha,ie),psi(1,1,alpha,ie),a,b,rofi,nr)
      norm = dsqrt(pr) ; psi(:,:,alpha,ie) = psi(:,:,alpha,ie)/norm

      enddo ! End of loop over energy (ie)
      enddo ! End of loop over alpha

C --- Orthonormalize Psi_1 and Psi_2
      do  ie = 1, 5
        pr = prodr(0,psi(1,1,1,ie),psi(1,1,2,ie),a,b,rofi,nr)
        psi(:,:,2,ie) = psi(:,:,2,ie) - pr * psi(:,:,1,ie)
C ...   Normalization of Psi_2 (missing in the previous versions)
        norm = dsqrt(prodr(0,psi(1,1,2,ie),psi(1,1,2,ie),a,b,rofi,nr))
        psi(:,:,2,ie) = psi(:,:,2,ie)/norm
      enddo

C --- Write wavefunctions to disk for testing purposes
c     if (l == 2) then
c       fil ='wfrel'
c       call awrit1('%a%i.dat',fil,11,0,imu)
c       print *,'filename:',fil, ':'
c       ifi = fopn(fil)
c       do  n = 1, nrmx
c         write(ifi,923)rofi(n),((psi(n,i,alp,3),i=1,4),alp=1,2)
c       enddo
c923    format(e12.4,8e15.6)
c       call fclr(fil,ifi)
c       print *,'Wrote (l,imu)=',l,imu
c       call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
c       read(*,'(a80)') outs
c       if (outs == 'q') call rx0('quit in rdeq')
c       if (imu == 2*(l+1)) stop
c     endif


C --- Save radial derivative at rmax ---
      if (imu == 1 .or. imu == 2*l+2) then
        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq
      else
        nug1(nr) = 1d0 + 2*z/rmax/csq + (enu1 - vm(nr) + up *bm(nr))/csq
      endif

      nug2(nr) = 1d0 + 2*z/rmax/csq + (enu2 - vm(nr) - upp*bm(nr))/csq

      gsmt(1,:) = (c*nug1(nr)*psi(nr,2,:,3) - xk /rmax*psi(nr,1,:,3))
      gsmt(2,:) = (c*nug2(nr)*psi(nr,4,:,3) - xk1/rmax*psi(nr,3,:,3))

      if (psi(nr,1,1,3)*psi(nr,3,2,3) < 0) then
        psi(:,:,2,:) = - psi(:,:,2,:)
      endif

C --- Save g and f at rmax ---
      gmt(1,:) = psi(nr,1,:,3)/rmax ; gmt(2,:) = psi(nr,3,:,3)/rmax
      fmt(1,:) = psi(nr,2,:,3)/rmax ; fmt(2,:) = psi(nr,4,:,3)/rmax

C --- Energy derivatives of rg,rf ---
      psd(:,:,:) = (-8 * (psi(:,:,:,2) - psi(:,:,:,4))
     .                  + psi(:,:,:,1) - psi(:,:,:,5))/12/de

C --- Orthogonalize phidot to phi ---
      do  alpha = 1, 2
C ...   Same alpha (already orthogonal; do to avoid roundoff error)
        pr = prodr(0,psi(1,1,alpha,3),psd(1,1,alpha),a,b,rofi,nr)
        psd(:,:,alpha) = psd(:,:,alpha) - pr * psi(:,:,alpha,3)
C ...   Different alphas
        pr = prodr(0,psi(1,1,3-alpha,3),psd(1,1,alpha),a,b,rofi,nr)
        psd(:,:,alpha) = psd(:,:,alpha) - pr * psi(:,:,3-alpha,3)
      enddo

C --- Save gdot and fdot at rmax ---
      gmtde(1,:) = psd(nr,1,:)/rmax; gmtde(2,:) = psd(nr,3,:)/rmax
      fmtde(1,:) = psd(nr,2,:)/rmax; fmtde(2,:) = psd(nr,4,:)/rmax

C --- Small parameter p=<gdot|gdot> by trapezium rule ---
      smalpa = 0
      do  n = 2, nrmx-1
        np = n + 1 ; r = rofi(n) ; r1 = rofi(np)
        w0 = (1 + xk *(xk+1) /((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
        w1 = (1 + xk *(xk+1) /((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
        w2 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r  - vm(n)) /c)*r)**2)
        w3 = (1 + xk1*(xk1+1)/((c + (eav + 2*z/r1 - vm(np))/c)*r1)**2)
        x1 = a*(r+b) ; x2 = a*(r1+b)
        do  a1 = 1, 2
          smalpa(a1,:) = smalpa(a1,:)
     .      + x1*(w0*psd(n ,1,a1)*psd(n ,1,:)+psd(n ,2,a1)*psd(n ,2,:))
     .      + x2*(w1*psd(np,1,a1)*psd(np,1,:)+psd(np,2,a1)*psd(np,2,:))
     .      + x1*(w2*psd(n ,3,a1)*psd(n ,3,:)+psd(n ,4,a1)*psd(n ,4,:))
     .      + x2*(w3*psd(np,3,a1)*psd(np,3,:)+psd(np,4,a1)*psd(np,4,:))
        enddo
        smalpa = smalpa/2
      enddo

      pprel(4,l,imu,:,:) = smalpa(:,:)

C --- Overlap integrals for the calculation of averages [Turek (6.141)]
C     Should be recorded for future use along with those for different mu's
C     (for now just print out)
C ... psi(n,kappa,alpha,ie)

c     write(6,900)l,mu
c900  format('Overlap parameters R for l=',i1,' , mu= ',F4.1)
c     do a1 = 1, 2
c       do a2 = 1, 2
c         do k1 = 1, 4
c           do k2 = 1, 4
c             if (mod(k1,2) /= mod(k2,2)) cycle ! don't need <g|f>
c             do ie = 1, 5
c              pr = prodr(1,psi(1,k1,a1,ie),psi(1,k2,a2,ie),a,b,rofi,nr)
c              if (mod(k1,2) == 1) then
c                write(6,901)ie,k1/2+1,a1,k2/2+1,a2,pr
c              else
c                write(6,902)ie,k1/2,a1,k2/2,a2,pr
c              endif
c             enddo
c           enddo
c         enddo
c       enddo
c     enddo
c901  format(1x,'ie=',i1,' Rgg(',i1,',',i1,';',i1,',',i1,')=',F11.8)
c902  format(1x,'ie=',i1,' Rff(',i1,',',i1,';',i1,',',i1,')=',F11.8)

      end

      subroutine fdpp(enu1,enu2,ksop,shft,gmt,fmt,gmtde,z,rmax,avw,l,lmx,imu,pprel)
C- Fully relativistic potential parameters
C------------------------------------------------------------------------------
Ci Inputs
Ci   enu1 :linearization energy, spin 1 (called enum1 in rdeq.f)
Ci   enu2 :linearization energy, spin 2 (called enum2 in rdeq.f)
Ci   ksop  :spin orbit coupling parameters, used to adjust enu1,enu2
Ci   shft  :constant shift of diagonal part of C parameter (diagonal part)
Ci   gmt   :2x2 normalized wave function g x r, at the MT surface
Ci   fmt   :2x2 normalized wave function f x r, at the MT surface
Ci   gmtde :gmtde :(2x2 matrix) gdot at rmax
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci         :used in the definition of potential parameters
Ci   l     :l quantum number
Ci   lmx   :dimensions ksop and pprel
Ci   imu   :l + mu + 3/2, mu = quantum number
Ci         :imu has range 1:2(l+1) and  mu has range (-l-1/2:l+1/2)
Co Outputs
Co   pprel :Relativistic potential parameters in kappa-mu representation
Co         :pprel = pprel(:,l=0:lmx,imu=1:2*(lmx+1),1:2,1:2)
Co         :Made here:
Co         :pprel(1,:,:,1:2,1:2) = C
Co         :pprel(2,:,:,1:2,1:2) = gamma
Co         :pprel(3,:,:,1:2,1:2) = D+, where (D+) D = delta
Co         :Not made here:
Co         :pprel(4,:,:,1:2,1:2) = small parameter p
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   18 Oct 13 (Belashchenko) cleanup
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer l,lmx,imu
      double precision gmt(2,2),fmt(2,2),gmtde(2,2),
     .                 mu, rmax, avw, z, enu1, enu2, shft
C ... Local parameters
      double precision pprel(5,0:lmx,2*(lmx+1),2,2)
      double precision gtmt(2,2), gmt1(2,2),
     .                 fg(2,2),
     .                 xk, xk1, kapa(2,2), D(2,2),
     .                 delta1(2,2),unit(2,2), q(2,2), q1(2,2),
     .                 gam(2,2), gam1(2,2), aa(2,2), aa1(2,2),
     .                 one(2,2), two(2,2), gmtdet(2,2), gamt(2,2),
     .                 yv(2,2), yv1(2,2), vm(2,2), enum(2,2),
     .                 gamma(2,2), gammat(2,2), delt(2,2),
     .                 deltt(2,2), cm(2,2), cm1(2,2),
     .                 savw, ksop(0:lmx,2,2,6),
     .                 xx(2,2),c
      real(8),parameter :: NULLR=-99999d0
C     double precision clebsh1(2,2), clebsh1t(2,2),clebsh(2,2), u,
C     double precision Ddot(2,2) Dtdot(2,2), gmtde1(2,2),fgde(2,2),delta(2,2),deltat(2,2),deltat1(2,2),fg1(2,2)

C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      if (z == 0) return

C --- Setup ---
      mu = dble(imu-l) - 1.5d0
      savw = (avw/rmax)**(2*l+1)

      if (imu == 1 .or. imu == 2*l+2) then
        xk1 = dble(-l) ; xk = xk1
      else
        xk1 = dble(-l) ; xk = dble(l+1)
      endif

      unit(1,1) = 1 ; unit(2,2) = 1 ; unit(1,2) = 0 ; unit(2,1) = 0

      enum = 0
      enum(1,1) = (enu1 + enu2)/2
      enum(2,2) = (enu1 + enu2)/2
      if (ksop(0,1,1,1) /= NULLR) then
        enum(1,1) = -(l+1) * ksop(l,1,1,1)/2 + (enu1 + enu2)/2
        enum(2,2) =   l    * ksop(l,1,1,1)/2 + (enu1 + enu2)/2
      endif
      if (imu == 1 .or. imu == 2*l+2) enum(1,1) = enum(2,2)

      kapa = 0 ; kapa(1,1) = xk ; kapa(2,2) = xk1

C --- Invert g and g_dot----------------------------
C     The g and f do not commute. In Phys Rev. B 43, 14414, 1991
C     it is  D = c*S*g^(-1)*f - k - I
C     In Sov. Phys. Solid State 31(8), 1285 1989 it is
C     D = c*S*f*g^(-1) - k - I.
C     Only the second formula gives symmetric D.
C     The order will depend on the order in g and f.
C     If alpha numerates columns then the second formula works,
C     If alpha numerates strokes then the first formula should be applied.
C----------------------------------------------------
      call dinv22(gmt,gmt1) ; fg = matmul(fmt,gmt1)
c     call dinv22(gmtde,gmtde1)
c     fgde = matmul(fmtde,gmtde1)

c     call mul22(gsmt,gmt1,fg1)

C --- Logarithmic derivative matrices-----------------
      D(:,:) = rmax*c*fg(:,:) - kapa(:,:)
c     D(:,:) = fg1(:,:) - unit(:,:)
c     Ddot(:,:) = rmax*c*fgde(:,:) - kapa(:,:)

C --- Symmetrize D (it should be symmetric, but for numerical reasons)
      D(1,2) = (D(1,2) + D(2,1))/2 ; D(2,1) = D(1,2)

C --- Make corrections for diagonal D, D_dot--------
c     if (mu < 0) then
c     D(1,1) = D(1,1) - srav1 + sr1
c     D(2,2) = D(2,2) - srav1 + sr1
c     Ddot(1,1) = Ddot(1,1) - sravdot1 + srdot1
c     Ddot(2,2) = Ddot(2,2) - sravdot1 + srdot1
c     else
c     D(1,1) = D(1,1) - srav2 + sr2
c     D(2,2) = D(2,2) - srav2 + sr2
c     Ddot(1,1) = Ddot(1,1) - sravdot2 + srdot2
c     Ddot(2,2) = Ddot(2,2) - sravdot2 + srdot2
c     endif

C------Make corrections for diagonal g, g_dot--------
c     if (mu < 0) then
c     gmt(1,1) = gmt(1,1) - gsrav1 + gsr1
c     gmt(2,2) = gmt(2,2) - gsrav1 + gsr1
c     gmtde(1,1) = gmtde(1,1) - gsravdot1 + gsrdot1
c     gmtde(2,2) = gmtde(2,2) - gsravdot1 + gsrdot1
c     else
c     gmt(1,1) = gmt(1,1) - gsrav2 + gsr2
c     gmt(2,2) = gmt(2,2) - gsrav2 + gsr2
c     gmtde(1,1) = gmtde(1,1) - gsravdot2 + gsrdot2
c     gmtde(2,2) = gmtde(2,2) - gsravdot2 + gsrdot2
c     endif

C --- Transpose gmt and Ddot -------------------------
      gtmt   = transpose(gmt) ; gmtdet = transpose(gmtde)
c     Dtdot = transpose(Ddot)

C --- Make unscreened potential parameters-------------
c     delta(:,:) = D(:,:) - Ddot(:,:)
c     deltat(:,:) = D(:,:) - Dtdot(:,:)
c     call dinv22(delta, delta1)
c     call dinv22(deltat, deltat1)
c     call mul22(gmtdet,gmt,one)
c     call mul22(gtmt,gmtde,two)
      call mul22(gmt,gmtdet,one) ! one = g x gdot^T
      call mul22(gmtde,gtmt,two) ! two = gdot x g^T

c     aa(:,:) = -(delta1(:,:)+deltat1(:,:))/2  ! Eq. (19) in Solovyev
      aa(:,:) = -rmax*(one(:,:) + two(:,:))/2  ! This symmetrization is superfluous
      call dinv22(aa,aa1)                      ! aa1 : A_L in Shick

      xx(:,:) = D(:,:) - l * unit(:,:)
      call dinv22(xx,delta1)                   !delta1 = (D-l)^-1
      q(:,:) = xx(:,:) + aa1(:,:)              !q = D + A - l
      call dinv22(q,q1)                        !q1 = (D + A -l)^-1
      q(:,:) = 2*(2*l+1) * savw * ((2*l+1)*q1(:,:)+unit(:,:))
                                               ! 2(2l+1) savw (D+A+l+1)(D+A-l)^-1
                                               ! q = gamma_L^-1 (inverse gamma from Shick)

      xx = matmul(aa,xx) + unit                ! xx = A^-1 * (D + A - l)
      call dinv22(xx,gam1)                     ! gam1 = (D + A - l)^-1 * A
      call mul22(gam1,gmt,xx)                  ! xx = (D + A - l)^-1 * A * g
      gam(:,:) = dsqrt(2d0*rmax*savw)*(2*l+1)*xx(:,:)
      gamt = transpose(gam)
      yv = delta1 + aa                         !yv = (D-l)^-1 + A^-1

      call dinv22(yv,yv1)                      !yv1 = ((D-l)^-1+A^-1)^-1
      call mul22(gtmt,yv1,yv)                  !yv = g^T * (..)
      call mul22(yv,gmt,vm)                    !vm = g^T * A [D + A - l]^-1 (D-l) * g

      vm(:,:) = rmax*vm(:,:) + enum(:,:)

C --- Potential parameters in the nearly orthogonal representation ---
      call dinv22(q,gamma)  ! this is gamma_L from Shick
      gammat = transpose(gamma) ! gamma should be symmetric anyway

cc       call mul22(gam,gamma,delt)
cc       call mul22(gammat,gamt,deltt)

      call mul22(gamma,gam,delt) ! delt = (..) * g
      call mul22(gamt,gamma,deltt) ! deltt = g^T * (..)

cc      call mul22(gam,gamma,cm)
cc     call mul22(cm,gamt,cm1)

cc       call mul22(delt,q,cm)
cc       call mul22(cm,deltt,cm1)

      call mul22(deltt,q,cm)  ! cm  = g^T * gam_L^-1
      call mul22(cm,delt,cm1) ! cm1 = g^T * (gam_L)^-1 * g

      cm(:,:) = vm(:,:) + cm1(:,:) + shft * unit(:,:)

      pprel(1,l,imu,:,:) = cm(:,:)
      pprel(2,l,imu,:,:) = gamma(:,:)
      pprel(3,l,imu,:,:) = delt(:,:)

C       open(190,file='cm.dat',status='unknown')
C       open(191,file='gamma.dat',status='unknown')
C       open(192,file='delt.dat',status='unknown')

c     do  i = 1, 2
c       do  j = 1, 2
C          write(190,*) cm(i,j)
C          write(191,*) gamma(i,j)
C          write(192,*) delt(i,j)
c       enddo
c     enddo

C --- Rotate cm---------------------------------------

cc       u = mu/(dble(l)+0.5d0)

cc       clebsh(1,1) = dsqrt(1d0+u)/dsqrt(2d0)
cc       clebsh(2,2) = clebsh(1,1)
cc       clebsh(2,1) = dsqrt(1d0-u)/dsqrt(2d0)
cc       clebsh(1,2) = -clebsh(2,1)

cc       call dinv22(clebsh,clebsh1)

cc      call mul22(clebsh1t,cm,cm1)
cc      call mul22(cm1,clebsh1,cm)

cc      if (dabs(mu) == (l+0.5d0)) then
cc      write(77,*) mu, xk-1, xk1-1
cc      write(77,*) '--------------------------------------------'
cc      write(77,*) cm(1,1), cm(1,2)
cc      write(77,*) cm(2,1), cm(2,2)
cc      write(77,*) '--------------------------------------------'
cc      endif


C-------This is a failed way to built the pot. pars-----------
cc       do  i = 1, 2
cc          do  j = 1, 2
cc             delta(i,j) = D(i,j) - Ddot(i,j)
cc          enddo
cc       enddo

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             deltat(i,j) = delta(j,i)
cc          enddo
cc       enddo


cc       call dinv22(delta, delta1)
cc       call dinv22(deltat, deltat1)

cc       call mul22(gmtdet,gmt,one)
cc       call mul22(gtmt,gmtde,two)

cc       do  i = 1, 2
cc        do  j = 1, 2
cc        aa(i,j) = -0.5d0*(delta1(i,j) + deltat1(i,j))
cc         aa(i,j) = -0.5d0*rmax*(one(i,j) + two(i,j))
cc        enddo
cc       enddo

cc       call dinv22(aa,aa1)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             delta(i,j) = D(i,j) - lmatrix(i,j)
cc           enddo
cc      enddo

cc       call dinv22(delta,delta1)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             q(i,j) = delta(i,j) + aa1(i,j)
cc          enddo
cc       enddo

cc      write(77,*) '--------------------------------------------'
cc      write(77,*) z, l, mu, xk, xk1
cc      write(77,*) '--------------------------------------------'
cc      write(77,*) delta(1,1), delta(1,2)
cc      write(77,*) delta(2,1), delta(2,2)
cc      write(77,*) '--------------------------------------------'
cc      write(77,*)  aa1(1,1), aa1(1,2)
cc      write(77,*)  aa1(2,1), aa1(2,2)

cc       call dinv22(q,q1)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             q(i,j) = 2*(2*l+1)*savw*((2*l+1)*q1(i,j)+unit(i,j))
cc          enddo
cc       enddo

cc       if (dabs(mu) == (l+0.5d0)) then
cc          write(77,*) l,q(1,1)
cc       endif

cc       call mul22(aa,delta,gam)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             gam(i,j) = gam(i,j) + unit(i,j)
cc          enddo
cc       enddo

cc       call dinv22(gam,gam1)
cc       call mul22(gam1,gmt,gam)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             gam(i,j) = dsqrt(2d0*rmax*savw)*(2*l+1)*gam(i,j)
cc          enddo
cc       enddo

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             gamt(i,j) = gam(j,i)
cc          enddo
cc       enddo

cc       do  i = 1, 2
cc          do  j = 1, 2
cc          yv(i,j) = delta1(i,j) + aa(i,j)
cc          enddo
cc       enddo

cc       call dinv22(yv,yv1)
cc       call mul22(gtmt,yv1,yv)
cc       call mul22(yv,gmt,vm)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc            vm(i,j) = rmax*vm(i,j) + enum(i,j)
cc          enddo
cc       enddo

C --- Potential parameters in the nearly orthonormalized representation
C------------------------------------------------------------

cc       call dinv22(q,gamma)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc             gammat(i,j) = gamma(j,i)
cc          enddo
cc       enddo

cc       call mul22(gam,gamma,delt)

cc       call mul22(gam,gamma,cm)
cc       call mul22(cm,gamt,cm1)

cc       do  i = 1, 2
cc          do  j = 1, 2
cc            cm(i,j) = vm(i,j) + cm1(i,j)
cc          enddo
cc       enddo

cc       open(190,file='cm.dat',status='unknown')
cc       open(191,file='gamma.dat',status='unknown')
cc       open(192,file='delt.dat',status='unknown')

cc       do  i = 1, 2
cc        do  j = 1, 2
cc           pprel(1,l,imu,i,j) = cm(i,j)
cc           write(190,*) cm(i,j)
cc           pprel(2,l,imu,i,j) = gamma(i,j)
cc           write(191,*) gamma(i,j)
cc           pprel(3,l,imu,i,j) = delt(i,j)
cc           write(192,*) delt(i,j)
cc        enddo
cc       enddo

C-------End of the failed way--------------------------------------

      end

      subroutine mul22(a,b,c)
C-    Multiplies 2x2 matrices a*b
C     a is from the left be carefull for non commuting matrices
      implicit none
      double precision a(2,2), b(2,2), c(2,2)

      !! ? replace with call to dmpy22

      c(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
      c(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
      c(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
      c(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

      end

C      subroutine product2(psi,ksi,a,b,rofi,nr,e,vm,z,xk,product)
CC-    Finds the product of two 4-vectors <psi|ksi> using the
CC     Trapezium rule
CC     Input parameters
C      integer nr
C      double precision psi(4,nr),ksi(4,nr),a,b,xk,rofi(nr),e,vm(nr),z
CC     Output
C      double precision product
CC     Local
C      integer n, m
CC     double precision cc
C
CC     cc = 274.072d0
C      product = 0d0
C
C      do  n = 2, nr-1
C        m = n + 1
C
CC         factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
CC     .                - vm(n))/cc)*rofi(n))**2)
CC         factor1 = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(m)
CC     .                - vm(m))/cc)*rofi(m))**2)
C
C        product = product + 0.5d0*a*((rofi(n)+b)*
C     .          ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n)) +
C     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n)) +
C     .         (rofi(m)+b)*
C     .          ((psi(1,m)*ksi(1,m) + psi(3,m)*ksi(3,m)) +
C     .            psi(2,m)*ksi(2,m) + psi(4,m)*ksi(4,m)))
C
C      enddo
C      end

      double precision function prodr(mode,w1,w2,a,b,rofi,nr)
      integer mode,nr
      double precision w1(nr,*),w2(nr,*),a,b,rofi(nr)
      integer i,imax,n
      double precision f1,f2

      if (mode == 0) then
        imax = 4
      elseif (mode == 1) then
        imax = 1
      else
        call rx('prodr: unknown mode')
      endif

      prodr = 0
      do i = 1, imax
        do  n = 2, nr-1
          f1 = a * (rofi(n)  + b)/2 ; f2 = a * (rofi(n+1) + b)/2
          prodr = prodr + f1*w1(n,i)*w2(n,i) + f2*w1(n+1,i)*w2(n+1,i)
        enddo
      enddo
      end

C      subroutine product4(psi,ksi,a,b,rofi,nr,e,vm,z,xk,product)
CC     Finds the product of two 4-vectors <psi|ksi> using the
CC     Simpson's rule
CC     Input parameters
C      integer nr
C      double precision psi(4,nr), ksi(4,nr), a, b, xk, rofi(nr), e,
C     .                 vm(nr), z
CC     Output
C      double precision product, product1, product2
CC     Local
C      integer n, m
C      double precision factor, factor1, cc
C
C      cc = 274.072d0
C      product = 0d0
C      product2 = 0d0
C      product1 = 0d0
C
C
C      do  n = 3, nr-2, 2
C
C          factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C
C         product = product + (rofi(n)+b)*
C     .            ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n))
C
C
C      enddo
C
C      product = (2.0d0*product)
C
C      do  n = 4, nr-1, 2
C
C         factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(n)
C     .                - vm(n))/cc)*rofi(n))**2)
C
C         product = product + (rofi(n)+b)*
C     .            ((psi(1,n)*ksi(1,n) + psi(3,n)*ksi(3,n))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,n)*ksi(4,n))
C
C         enddo
C
C
C      factor = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(2)
C     .                - vm(2))/cc)*rofi(2))**2)
C
C      factor1 = (1d0 + xk*(xk+1d0)/((cc + (e + 2d0*z/rofi(nr)
C     .                - vm(nr))/cc)*rofi(nr))**2)
C
C      product1 = (rofi(2)+b)*
C     .            ((psi(1,2)*ksi(1,2) + psi(3,2)*ksi(3,2))*factor +
C     .            psi(2,n)*ksi(2,n) + psi(4,2)*ksi(4,2))
C
C      product2 = (rofi(nr)+b)*
C     .            ((psi(1,nr)*ksi(1,nr) + psi(3,nr)*ksi(3,nr))*factor1 +
C     .            psi(2,nr)*ksi(2,nr) + psi(4,nr)*ksi(4,nr))
C
C
C      product = (2.0d0*product + product1 + product2)*a/3.0d0
C
C      return
C      end

      subroutine clebsh(l,imu,a,a1)
      implicit none
      double precision a(2,2),a1(2,2),cl(2,2),mu,u
C     double precision cl1(2,2)
      integer ms1,ms2,l,imu


      mu= dble(imu - l) - 1.5d0
      u = mu/(dble(l)+0.5d0)

      cl(1,1) = dsqrt(0.5d0*(1d0+u))
      cl(2,2) = cl(1,1)
      cl(2,1) = -dsqrt(0.5d0*(1d0-u))
      cl(1,2) = -cl(2,1)

C         call dinv22(cl,cl1)

      do  ms1 = 1, 2
        do  ms2 = 1, 2
          a1(ms1,ms2) =
     .      cl(1,ms1)*a(1,1)*cl(1,ms2) +
     .      cl(1,ms1)*a(1,2)*cl(2,ms2) +
     .      cl(2,ms1)*a(2,1)*cl(1,ms2) +
     .      cl(2,ms1)*a(2,2)*cl(2,ms2)
        enddo
      enddo

C       --- Test rotate back to lms rep'sn from kappa-mu rep'sn  ---
C          do  ms1 = 1, 2
C             do  ms2 = 1, 2
C                a(ms1,ms2) =
C     .            cl1(1,ms1)*a1(1,1)*cl1(1,ms2) +
C     .            cl1(1,ms1)*a1(1,2)*cl1(2,ms2) +
C     .            cl1(2,ms1)*a1(2,1)*cl1(1,ms2) +
C     .            cl1(2,ms1)*a1(2,2)*cl1(2,ms2)
C
C              else
C                reslms(m1,m2,ms1,ms2) = 0d0
C              endif
C
C            enddo
C            enddo

      end
