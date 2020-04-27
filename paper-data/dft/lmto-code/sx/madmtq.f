      subroutine madmtq(mode,q,nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg
     .  ,cmad,v1,v2)
C- Madelung matrix for nonzero q
C ----------------------------------------------------------------
Ci Inputs
Ci   mode  0, generate cmad(q)
Ci         1, aproximate cmad to leading order in q, with the singular
Ci            term (v2/q**2) subracted (NB : imaginary part untouched)
Ci            and make coefficients v1,v2
Ci   nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg
Co Outputs
Co   cmad     Madelung matrix (mode=0) or an estimate of
Co            the same, with the 1/q**2 terms subtracted (mode=1)
Co   v1       (mode=1) coff to q.tau term in cmad; see Remarks
Co   v2       (mode=1) coff to 1/q**2 term in cmad; see Remarks
Cr Remarks
Cr   Madelung matrix has the limiting behavior as q->0
Cr     M(q->0) = M0 + v2/q**2 - 1/2 v2 [q.tau/|q|/2/pi]**2
Cr               where M0 = result of strxq0 for q=0
Cr             = M0 + v2/q**2 +/- 1/2 v1 [q.tau/|q|/2/pi]
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,nkd,nkg
      double precision awald,alat,vol,q(3)
      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg),
     .  cmad(nbas,nbas,2),v1(nbas,nbas),v2
C ... Local parameters
      integer ibas,jbas,iprint,i,lgunit
      double precision tau(3),tau0(3),cxx(2),q0(3),pi,qq,qr,ddot

      if (mode == 1) then
        call dpzero(q0,3)
        pi = 4.d0*datan(1.d0)
        v2 = alat**2/vol/pi
      else
        call dcopy(3,q,1,q0,1)
      endif

C --- Generate Madelung matrix ---
      do  10  ibas = 1, nbas
      do  10  jbas = ibas, nbas
        do  15  i = 1, 3
   15   tau0(i) = bas(i,ibas)-bas(i,jbas)
        call shortn(tau0,tau,dlat,nkd)
        call strxq0(tau,tau0,awald,alat,vol,glat,nkg,dlat,nkd,cxx,q0)
        cmad(ibas,jbas,1) = cxx(1)
        cmad(jbas,ibas,1) = cxx(1)
        if (mode == 0) then
          cmad(ibas,jbas,2) =  cxx(2)
          cmad(jbas,ibas,2) = -cxx(2)
        else
          qq = ddot(3,q,1,q,1)
          qr = ddot(3,tau0,1,q,1)/dsqrt(qq)*2d0*pi
          v1(ibas,jbas) =  qr*v2
          v1(jbas,ibas) = -qr*v2
          cmad(ibas,jbas,1) = cmad(ibas,jbas,1) - 0.5d0*qr**2*v2
          cmad(jbas,ibas,1) = cmad(jbas,ibas,1) - 0.5d0*qr**2*v2
        endif
   10 continue
      if (mode /= 0) v2 = 2d0*v2


C --- Printout ---
      if (iprint() < 60) return
      call awrit2('%N MADMTQ: Madelung matrix%?;n; V0;;, q =%3:1,3;3d'
     .  ,' ',80,lgunit(1),mode,q)
      do  20  ibas = 1, nbas
   20 print 110, (cmad(ibas,jbas,1),jbas=1,nbas)
      if (mode == 0) then
        print *, 'Imaginary part'
        do  22  ibas = 1, nbas
   22   print 110, (cmad(ibas,jbas,2),jbas=1,nbas)
      else
        call awrit1(' v2=%,4;4d  v1:',' ',80,lgunit(1),v2)
        do  24  ibas = 1, nbas
   24   print 110, (v1(ibas,jbas),jbas=1,nbas)
      endif
  110 format(9f8.3)
      end

      subroutine strxq0(tau,tau0,awald,alat,vol,glat,nkg,dlat,nkd,dl,q)
C- Widget to make structure constant DL for L=0,E=0,K=0.
C ----------------------------------------------------------------
Ci Inputs
Ci   TAU,awald,ALAT,
Ci   VOL:     cell volume
Ci   glat,nkg:reciprocal lattice vectors and number
Ci   DLAT,NKD:real (direct) space lattice vectors and number
Co Outputs
Co   dl: 'potential' Phi at tau (see remarks below)
Cr Remarks
Cr   dl is in inverse atomic units.
Cr
Cr  The 1/r structure constant is calculated by means of a Madelung sum
Cr  of a set of discrete point charges 1/r compensated by a background
Cr  of uniform density that makes total charge neutral.  The sum of
Cr  'potentials' 1/r at the position tau is returned as DL.  The
Cr  potential from a basis of more than one atom is obtained by
Cr  superposing the potential from each sublattice.
Cr
Cr  The potential is broken into a "damped" and a "smooth" part as
Cr  (1) h(r) = 1/r = erfc(a*r)/r + erf(a*r)/r = h^d(r) + h^s(r)
Cr  which is independent of the Ewald parameter a, and the sum
Cr  (2) Phi(r) =
Cr  is calculated by adding contributions from h^d and h^s separately;
Cr  by subtracting the uniform background total potential is finite.
Cr  In the case r -> 0, The onsite term h(r) is subracted making the
Cr  potential finite at the origin.
Cr  h^s is calculated in Fourier space as
Cr  (3)
Cr     4
Cr  The term G=0 is infinite in consequence of the infinite range of
Cr  the 1/r potential, but is canceled by a uniform background, up to
Cr  a constant of integration.  The constant is chosen to make Phi
Cr  independent of a, and also to make Phi -> 0 as the density of
Cr  lattice points goes to zero.  Equation (2) becomes:
Cr  (4) Phi(r) =
Cr
Cr
Cr  and in the special case r->0 it is noted that
Cr  (5)
Cr  and the contribution h^s(0) = 2*a/
Cr
Cr  First lattice vector must be the zero vector.
C ----------------------------------------------------------------
C      implicit none
C Passed parameters
      integer nkg,nkd
      double precision tau(3),glat(3,nkg),dlat(3,nkd),q(3),tau0(3)
      double precision awald,alat,vol
      double complex dl
c local parameters
      integer ir,ir1
      double precision pi,tpi,gamma,tpiba2,r1,r2,derfc,d1mach,sp
      external derfc,d1mach

      pi = 4.d0*datan(1.d0)
      tpi = 2.d0*pi
      gamma = 0.25d0/awald**2
      tpiba2 = (tpi/alat)**2
C --- Reciprocal space sums ---
      if (q(1)**2 + q(2)**2 + q(3)**2 > d1mach(3)) then
        ig1 = 1
        dl=0d0
      else
        ig1 = 2
        dl = -gamma
      endif
      do  26  ir = ig1,nkg
        r2 = tpiba2*((q(1)+glat(1,ir))**2 + (q(2)+glat(2,ir))**2 +
     .               (q(3)+glat(3,ir))**2)
        sp =tpi*((q(1)+glat(1,ir))*tau0(1)+ (q(2)+glat(2,ir))*tau0(2)+
     .           (q(3)+glat(3,ir))*tau0(3))
        dl = dl + dexp(-gamma*r2)/r2*dcmplx(cos(sp),sin(sp))
   26 continue
      dl = dl*4.d0*pi/vol

C --- Real space sums ---
      if (tau(1)**2 + tau(2)**2 + tau(3)**2 > d1mach(3)) then
        ir1 = 1
      else
        ir1 = 2
        dl = dl - 2.d0*awald/dsqrt(pi)
      endif
      do  20  ir = nkd, ir1, -1
        sp=tpi*(q(1)*(dlat(1,ir)+tau0(1)-tau(1))+
     .          q(2)*(dlat(2,ir)+tau0(2)-tau(2))+
     .          q(3)*(dlat(3,ir)+tau0(3)-tau(3)))
        r1 = alat*dsqrt((tau(1)-dlat(1,ir))**2 +
     .                  (tau(2)-dlat(2,ir))**2 +
     .                  (tau(3)-dlat(3,ir))**2)
        dl = dl + derfc(awald*r1)/r1*dcmplx(cos(sp),sin(sp))

c|      write(*,'(a,3f10.3)') ' r1,sp:',r1/alat,sp/tpi

   20 continue
      return
      end
