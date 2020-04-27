      subroutine enutod(lrelon,nl,l,ic,rmax,avw,pp,amom,amomr,vdel,idmod,
     .  pmin,pnu,qnu,qnur,eb)
C- Translate band-center energy and shift moments for one l
C ----------------------------------------------------------------
Ci Inputs
Ci   lrelon:T => fully relativistic mode
Ci   l     :quantum number with moments
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   amom: moments corresponding to potential parameters pp
Ci   vdel: difference in origins of energy scales by which the moments
Ci         and enu's are defined:  moments and enu defined with a
Ci         common origin => vdel is 0
Ci   eb:   (idmod=2 only): prescribed shift in enu
Ci   idmod:0 => pnu, qnu are shifted to center of gravity of band
Ci         1 => pnu is frozen; qnu shifted only by vdel.
Ci         2 => pnu, qnu to make change in enu=eb+vdel
Ci         3 => pnu, qnu floated to enu = C
Co Outputs
Co   pnu   :may be shifted depending on idmod
Co   qnu   :shifted to make moments relative to pnu;
Co   qnur  :also shifted if lrelon is T
Co   eb    :energy shift corresponding to shift in qnu and pnu.
Cr Remarks
Cr   Principal quantum number = integer part of PNU is retained.
Cr   amom and qnu can point to the same address space.
Cr   Note: should be checked for vdel /= 0
Cr
Cr   Relation between E and D as follows:
Cr   omeg2 - omeg1 = -S Phi1 Phi2 (D2 - D1), in general,
Cr   so ... for omeg = e - enu, find D(omeg) as follows:
Cr   omeg - omegp = -S Phip Phi(omeg) (D(omeg) - Dp)
Cr   omeg - omegm = -S Phim Phi(omeg) (D(omeg) - Dm)
Cr   Use these to eliminate Phi(omeg) and recall that
Cr   Phim = (2 delta/S sdivw)^1/2 and Phip = Phim sdivw/(2(2l+1) Gamma)
Cr   and  omegm = C - enu and omegp = omegm - delta/Gamma
Cr
Cr   Transformation of the moments:
Cr   The zeroth "moment" m0, i.e. int. phi^2 dE, is related to the
Cr   input moment q0 by
Cr     (1) q0 = m0 + p m2.
Cr   with q0 equal to the the total charge inside the sphere.
Cr   The other moments q1 and q2 are identical to m1 and m2.
Cr
Cr   When shifting by e = change in enu, the wave functions
Cr   phi and phidot are rotated to phit and phitdot as:
Cr     (2) phi    = (phit - e*phitdot)/n,
Cr         phidot = p*e*phit + phitdot)/n
Cr   where
Cr     (3) n^2 = 1 + p*e**2
Cr   This transformation preserves normalization of phi and phidot
Cr   and orthogonality of phi with phidot.  The charge density is
Cr     (4) n(r) := m0*phi**2 + 2*m1*phi*phid + m2*phid**2;
Cr   and the charge density in the shifted representation is
Cr   obtained from the linear transformation (2).  The coefficients
Cr   q0,q1 and q2 corresponding to phit and phitdot are obtained
Cr   from the coefficients q0,q1,q2 related to phi and phidot as:
Cr
Cr     (5) q0t = q0
Cr                    2                                    2
Cr         q1t =  (- e *p*q1 + 2*e*p*q2 - e*q0 + q1)/(1 + e *p)
Cr                    2         2                         2
Cr         q2t =  (- e *p*q2 + e *q0 - 2*e*q1 + q2)/(1 + e *p)
Cr
Cr   When corrections order p are neglected:
Cr
Cr   q0t = q0,  q1t = q1 - e q0,  q2t = e^2 q0 - 2 e q1 + q2
Cu Updates
Cu   21 Apr 15 (Vishina) extended to relativistic pp's
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer idmod,l,nl,ic
      logical lrelon
      double precision vdel,rmax,pmin
      double precision amom(3),pp(6),pnu,qnu(3),avw,eb
      double precision qnur(4,0:nl-1,2*nl,2,2),amomr(4,0:nl-1,2*nl,2,2)
C ... Local parameters
      logical lpmin
      integer ipqn,imu,i1,i2,ipass
      double precision omeg,omegm,omegp,phim,phip,aplus,amins,dnu,sdivw
      real(8), parameter :: pi=4*datan(1d0)
      double precision q0,q1,q2,d1mach,pfree,dmin,pfuni
      double precision q0r(2*nl,2,2),q1r(2*nl,2,2),q2r(2*nl,2,2),ebr(2*nl)
C     integer procid,master,mpipid

C      master = 0
C      procid = mpipid(1)

      ipqn = pnu
      if (lrelon) ebr = 0
      if (mod(idmod,10) /= 2 .or. dabs(amom(1)) <= d1mach(3)) eb = 0

C Check if no charge in this channel
      if (dabs(amom(1)) <= d1mach(3)) then
        call info5(10,0,0,' enutod (warning) ic=%i l=%i q0=%,4;4g < 0 ... no shift',ic,l,amom(1),4,5)
        return
      endif

C Shift enu to c.g. if idmod is zero; otherwise shift moments to enu
      if (mod(idmod,10) == 0) then
        eb = amom(2)/amom(1)
      elseif (mod(idmod,10) == 3) then
        eb = pp(2)-pp(1)
      else
        eb = eb+vdel
      endif

C Skip estimate of pnu if freeze pnu and also already estimate
      if (pnu-ipqn > .001 .and. mod(idmod,10) == 1) goto 20

      pfree = 0.5d0 - datan(dble(l))/pi
      if (pmin == 1) pmin = pfree
      if (pmin > 0) pmin = ipqn + pmin

      ipass = 0
   10 continue
      ipass = ipass+1

C     omega(E) = E(P)-enu = C-enu + [(P^-1-gamma)/delta]^-1
C     omegm = omega(D=-l-1) = omega(P=0) = C-enu
C     omegp = omega(D=l) = omega(P^-1=0) = C-enu-delta/gamma
      sdivw = (rmax/avw)**(l+l+1)
      omeg = eb - vdel
      omegm = pp(2) - pp(1)
      omegp = omegm - pp(3)**2/pp(5)
      phim = dsqrt(2/rmax) * pp(3)/dsqrt(sdivw)
      phip = phim / (2*(2*l+1)*pp(5)/sdivw)
      aplus = (omeg-omegp)/phip
      amins = (omeg-omegm)/phim
      dnu = (amins*l + aplus*(l+1)) / (amins-aplus)
      pnu = .5d0 - datan(dnu)/pi + ipqn
C     Case pnu < pmin:  find eb corresponding to pmin
C     To second order, E(P)=enu+omega(D)
      lpmin = .false.
      if (pnu < pmin) then
        lpmin = .true.
        dmin = -tan(pi*(pmin-ipqn-.5d0))
C       Inverse potential function
        pfuni = sdivw/(2*(2*l+1))*(dmin-l)/(dmin+l+1)
C       P^-1(E) = delta/(E-C)+gamma ->
C       E-enu = C-enu + [delta/(P^-1-gamma)]
        eb = pp(2)-pp(1) + pp(3)*pp(3)/(pfuni-pp(5))
        pnu = pmin
      endif

C --- Transform moments ---
   20 continue

C ... Find enu shift for relativistic moments
      if (lrelon) then
        do  imu = 1, 2*l+2
          if (dabs(amomr(1,l,imu,1,1)+amomr(1,l,imu,2,2)) >= d1mach(3)) then
            ebr(imu) = (amomr(2,l,imu,1,1)+amomr(2,l,imu,2,2))/
     .                 (amomr(1,l,imu,1,1)+amomr(1,l,imu,2,2))
          else
            ebr(imu) = 0
          endif
        enddo
      endif

      if (amom(3) < 0) then
        call info5(10,0,0,' enutod (warning) ic=%i l=%i q2=%,4;4g < 0 ... reset to 0',ic,l,amom(3),4,5)
        amom(3) = 0
      endif
      q0 = amom(1)
      q1 = amom(2) - eb*amom(1)
      q2 = amom(3) - eb*(amom(2) + q1)
      if (lrelon) then
        q1r = -9999
        do imu = 1, 2*l+2
          do i1 = 1,2
          do i2 = 1,2
           q0r(imu,i1,i2) = amomr(1,l,imu,i1,i2)
           q1r(imu,i1,i2) = amomr(2,l,imu,i1,i2) - ebr(imu)*amomr(1,l,imu,i1,i2)
           q2r(imu,i1,i2) = amomr(3,l,imu,i1,i2) - ebr(imu)*(amomr(2,l,imu,i1,i2)+q1r(imu,i1,i2))
          enddo
          enddo
        enddo
      endif

C ... Handle special case quadratic approximation to shifts makes q2<0
C     Note: this will never happen if eb was zero in the first place
      if (q2 < 0 .and. lpmin) then ! Already at boundary; best we can do is set q2=0
        call info5(10,0,0,' ENUTOD ic=%i l=%i P=pmin=%,4;4d but q2=%,4;4g.  Use EB=%,4;4d and set q2=0',
     .    ic,l,pmin,q2,eb)
        q2 = 0
      elseif (q2 < 0) then
        omegp = -(sqrt(amom(2)**2-amom(1)*amom(3))-amom(2))/amom(1)
        omegm = (sqrt(amom(2)**2-amom(1)*amom(3))+amom(2))/amom(1)
C       pass 1 : shift eb to closest branch that makes q2=0
        if (ipass <= 2) then
C         if (abs(eb-omegp) < abs(eb-omegm)) then
          if (abs(omegp) < abs(omegm)) then
            omeg = omegp*.999d0
C           print *, 'branch 1',omegp,omegm
          else
            omeg = omegm*.999d0
C           print *, 'branch 2',omegp,omegm
          endif
          if (ipass == 2) omeg = omeg/2
        else ! Give up ... revert to no shift
          omeg = 0
        endif
        if (amom(3) == 0) omeg = 0  ! unshifted moment already <=0 ... be conservative
        aplus = amom(2) - omeg*amom(1)
        amins = amom(3) - omeg*(amom(2) + aplus)
        call info8(10,0,0,' ENUTOD ic=%i l=%i eb=%,4;4d q2=%,4;4g.  Try EB=%,4;4d for q2=%,4;4g',
     .    ic,l,eb,q2,omeg,amins,7,8)
C       print *, 'ipass',ipass,lpmin,amom(3),q2
        eb = omeg
        goto 10
      endif
      qnu(1) = q0
      qnu(2) = q1
      qnu(3) = q2

      if (lrelon) then
        do  imu = 1, 2*l+2
          qnur(1,l,imu,:,:) = q0r(imu,:,:)
          qnur(2,l,imu,:,:) = q1r(imu,:,:)
          qnur(3,l,imu,:,:) = q2r(imu,:,:)
          qnur(4,l,imu,:,:) = qnur(4,l,imu,:,:) + ebr(imu)
        enddo
      endif

c     norm = 1 + eb**2*pp(4)
c     n2 =   1 - eb**2*pp(4)
c     qnu(1) = q0
c     qnu(2) = (2*eb*pp(4)*q2 - eb*q0 + q1*n2)/norm
c     qnu(3) = (eb**2*q0 - 2*eb*q1 + q2*n2)/norm
      end

      subroutine qrel2z12(nl,nc,qnur)
C- Zero out coefficient to relativistic phidot(1)*phidot(2)
      implicit none
C ... Passed parameters
      integer nl,nc
      double precision qnur(0:3,0:nl-1,2*nl,2,2,nc)
C ... Local parameters
      integer ic

      do  ic = 1, nc
        qnur(2,0:nl-1,1:2*nl,1,2,ic) = 0
        qnur(2,0:nl-1,1:2*nl,2,1,ic) = 0
      enddo

      end

      subroutine qtotrel(nl,nc,qnur,qtot)
C- Return total charge from relativistic moments
      implicit none
C ... Passed parameters
      integer nl,nc
      double precision qnur(4,0:nl-1,2*nl,2,2,nc),qtot(2,nc)
C ... Local parameters
      integer ic,l,imu,i1,i2
      double precision mu,Mrel

      do  ic = 1, nc

C       Charge
        qtot(1,ic) = sum(qnur(1,:,:,1,1,ic))+sum(qnur(1,:,:,2,2,ic))

C       Magnetic moment
        Mrel = 0
        do  l = 0, nl-1
        do  imu = 1, 2*(l+1)
          mu  = imu - l - 1.5d0
          if (imu == 1 .or. imu == 2*l+2) then
            Mrel = Mrel -2*mu*qnur(1,l,imu,2,2,ic)
          else
            do i1 = 1,2
              do i2 = 1,2
                Mrel = Mrel - 2*mu*qnur(1,l,imu,i1,i2,ic)
              enddo
            enddo
          endif
        enddo
        enddo
        qtot(2,ic) = Mrel
      enddo

      end
