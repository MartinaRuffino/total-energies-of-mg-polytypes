      subroutine rhocor(isw,z,lmax,nsp,konfig,a,b,nr,rofi,v,g,kcor,lcor,
     .                  qcor,tol,ec,sumec,sumtc,rho,gcore,ipr)
C- Generates the (spherical) charge density from the core states
C ----------------------------------------------------------------------
Ci Inputs:
Ci   isw   :1s digit
Ci         :1 return the core wave functions in gcore
Ci         :10s digit
Ci         :0 Use v for potential
Ci         :1 Use for potential average of spin-up, spin-down
Ci         :  Then spin-dependence of core eigenenergy originates
Ci         :  from potential due to valence electron only.
Ci         :100s digit
Ci         :0 Treat core scalar relativistically
Ci         :1 Calculate core wave functions using Dirac equation
Ci   z     :nuclear charge
Ci   lmax  :maximum l
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci   konfig:core configuration. Core orbitals are specified by:
Ci         :  1, 2, ..., konf(0)-1 for s
Ci         :  2, 3, ..., konf(1)-1 for p
Ci         :  3, 4, ..., konf(2)-1 for d, and so on.
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (electronic contribution)
Ci   g     :work array holding normalized wave function times r
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   tol   :wave function tolerance
Ci   ipr   :0 no printout
Ci         :1 summary printout
Ci         :2 detailed printout
Cio Inputs/Outputs:
Cio  ec    :guessed core eigenvalues (input)
Cio         core eigenvalues (output)
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co         :the spherical charge density of the core states is added
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy
Co   gcore :(isw=0) not used
Co         :(isw=1) core wave functions
Cl Local variables
Cl   deg   :orbital occupation number
Cb Bugs
Cb   calls info ... should not do this since may be called in parallel!
Cr Remarks:
Cu Updates
Cu   20 Aug 16 First cut at fully relativistic core
Cu   28 Jan 13 core can be calculated with spin-averaged potential
Cu   19 Apr 01 core wave functions may be saved in gcore
C  ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer isw,ipr,lmax,konfig(0:lmax),nr,nsp,kcor,lcor
      double precision a,b,z,g(nr,2),rho(nr,nsp),rofi(nr,2),tol,qcor(2),
     .  sumec(nsp),sumtc(nsp),gcore(nr,2,*),ec(*)
      real(8), target  :: v(nr,nsp)
C ... Dynamically allocated arrays
      real(8), pointer  :: vloc(:,:)
      real(8), allocatable  :: rhoD(:,:)
C     real(8), allocatable  :: psi(:,:,:)
C ... Local parameters
      character*1 pqn(9),ang(6), outs*128
      integer icore,ncore,intopt,isp,isw0,isw1,isw2,kc,konf,l,il,nodes,
     .  nre,stdo,kpr,nreD(0:lmax,0:lmax)
      double precision deg,dlml,e1,e2,ecor0,ecore,qcore,rhorim,rmax,
     .  rorim,slo,sum,tcore,val,vrho,kappa2,tc(2*(lmax+1)**2)
C     For Dirac equation
      double precision c,E,rorimD,sumecD,sumtcD,xv(4),
     .  Ekmu(2*(2*lmax+1),0:lmax,0:lmax),tcorD(0:lmax,0:lmax)
      procedure(integer) :: nglob
      procedure(real(8)) :: dsum,ddot

C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

C ... Heap
      data pqn /'1','2','3','4','5','6','7','8','9'/
      data ang /'s','p','d','f','g','h'/

C     print *, '!! setting v=0'; v=0

      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      isw0 = mod(isw,10)
      isw1 = mod(isw/10,10)
      isw2 = mod(isw/100,10)
      sumecD = 0; sumtcD = 0; rorimD = 0

      kpr = ipr+2*0
      if (isw2 /= 0) kpr = -kpr  ! print both full and scalar Dirac data

C     Generate vloc = potential with 2 spin channels and B=0
C            (spin-average )               (spin-split)
      if (isw1 /= 0.and.nsp == 2 .or. isw2 /= 0.and.nsp == 1) then
        allocate(vloc(nr,2))
        call dpscop(v,vloc,nr,1,1,0.5d0)
        call daxpy(nr,0.5d0,v(1,nsp),1,vloc,1)
        call dcopy(nr,vloc,1,vloc(1,2),1)
C     Otherwise use v
      else
        vloc => v
      endif
      if (isw2 /= 0) then
        allocate(rhoD(nr,2))
        call dpzero(rhoD,nr*2)
      endif

      rmax  = rofi(nr,1)
      b     = rmax/(dexp(a*nr-a)-1d0)
C     call radwgt(intopt,rofi(nr,1),a,nr,rofi(1,2))
      e1    = -2.5d0*z*z - 5d0
      e2    = 20d0
      icore = 0
      do  isp = 1, nsp
        sumec(isp) = 0d0
        sumtc(isp) = 0d0
        qcore = 0d0
        rhorim = 0d0
        if (kpr >= 2) write(stdo,1)
        if (kpr >= 2.and.nsp == 2) write(stdo,'('' spin'',i2,'':'')')isp
    1   format(/' state  chg          ecor0',10x,'ecore',10x,'tcore',
     .     4x,'nre',2x,'rho(rmax)')
        do  l = 0, lmax
          do  konf = l, konfig(l)-2
            deg = (2*(2*l+1))/nsp
            if (konf+1 == kcor .and. l == lcor) then
              deg = deg + (qcor(1)+(3-2*isp)*qcor(2))/nsp
            endif
            icore = icore+1
            nodes = konf - l
            ecor0 = ec(icore)
            val = 1d-30
            slo = -val
            call rseq(e1,e2,ecor0,tol,z,l,nodes,val,slo,vloc(1,isp),
     .        g,sum,a,b,rofi,nr,nre,kc)
C            print *, 'l, nodes',l,nodes
C            call prrmsh('g(SR)',rofi,g,nr,nr,2)
            ecore = ecor0

C       ... Correct core energy by using hankel bc's
            kappa2 = ecor0 - vloc(nr,isp) + 2*z/rmax
            if (nre == nr .and. kappa2 < 0d0) then
              dlml = -1d0 - dsqrt(-kappa2)*rmax
                do  il = 1, l
                  dlml = -kappa2*rmax*rmax/dlml - (2*il+1)
                enddo
              slo = val*(dlml+l+1)/rmax
              call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,vloc(1,isp),
     .          g,sum,a,b,rofi,nr,nre,kc)
            endif
            ec(icore) = ecore
            if (isw0 == 1) call dcopy(2*nr,g,1,gcore(1,1,icore),1)

C       ... Scalar Dirac: add to rho, make integral v*rho
C           ? why not use vloc(1,isp)
            call xyrhsr(ecore,l,z,nr,nre,g,rofi,rofi(1,2),
     .        vloc(1:nr,isp),rho(1,isp),deg,vrho,rorim)
            tcore = ecore - vrho
            qcore = qcore + deg
            rhorim = rhorim + rorim
            tc(icore) = tcore
            sumec(isp) = sumec(isp) + deg*ecore
            sumtc(isp) = sumtc(isp) + deg*tcore
            if (kpr >= 2) write (stdo,2) pqn(konf+1),ang(l+1),deg,
     .        ecor0,ecore,tcore,nre,rorim
    2       format(1x,2a1,f8.2,3f15.6,i7,f9.5)

C       ... Full Dirac equation: add to relativistic rho, make v*rho
            if (isw2 /= 0 .and. isp == 1) then
              if (kcor /= 0)
     .          call rx('rhocor: Dirac core not ready for partial core occupation')
              E = ecore
C             allocate(psi(nr,4,2*(2*l+1)))
              call rdeqcore(10+isw2,E,z,vloc,rofi,nr,nsp,nodes,kc,a,b,l,
     .          Ekmu(1,l,konf),rhoD,tcorD(l,konf),nreD(l,konf),rorimD)
              sumecD = sumecD + dsum(2*(2*l+1),Ekmu(1,l,konf),1)
              sumtcD = sumtcD + tcorD(l,konf)
C             deallocate(psi)
C             print *, 'XX',tcore*deg - tcorD(l,konf), sumtcD-sumtc(1)
            endif

          enddo                 ! konf
        enddo                   ! l
        if (kpr > 0) write(stdo,3) qcore,sumec(isp),sumtc(isp),rhorim
    3   format(' sum q=',f5.2,'  sum ec=',f15.8,'  sum tc=',f15.8,'  rho(rmax)',f8.5)
      enddo ! spin
      ncore = icore/2; if (nsp == 1) ncore = 0

      if (isw2 /= 0 .and. kpr <= -2) then
        write(stdo,4)
    4   format(/' Dirac core levels:'/' nl  chg',
     .    4x,'<ecore(S)>',5x,'<ecore(D)>',
     .    5x,'<Tcore(S)>',5x,'<Tcore(D)>',
     .    3x,'nre')
C    .    3x,'nre',2x,'rho(rmax)')

        icore = 0
        do  l = 0, lmax
          do  konf = l, konfig(l)-2
            icore = icore+1
            write (outs,2) pqn(konf+1),ang(l+1)
            call info8(1,0,0,outs//'%a%,4i%;15,6D%;15,6D%;15,6D%;15,6D%,6i',(2*(2*l+1)),
     .        (ec(icore)+ec(icore+ncore))/2,
     .        dsum(2*(2*l+1),Ekmu(1,l,konf),1)/(2*(2*l+1)),
     .        (tc(icore)+tc(icore+ncore))/2,
     .        tcorD(l,konf)/(2*(2*l+1)),nreD(l,konf),7,8)
            call info2(1,0,0,' ec(mu)%n:;14,6D',2*(2*l+1),ekmu(1,l,konf))
          enddo
        enddo
      endif
      if (isw2 /= 0 .and. kpr <= -1) then
        call info5(1,1,0,' qcore(SR) %;9F  qcore(FR)  %;9F  rho(rmax) %;8,5D',
     .    ddot(nr,rho,1,rofi(1,2),1)+(nsp-1)*ddot(nr,rho(1,nsp),1,rofi(1,2),1),
     .    ddot(nr,rhoD,1,rofi(1,2),1)+ddot(nr,rhoD(1,2),1,rofi(1,2),1),
     .    rhoD(nr,1)+rhoD(nr,2),4,5)
        xv(1) = sumec(1)+(nsp-1)*sumec(nsp)
        xv(2) = sumtc(1)+(nsp-1)*sumtc(nsp)
        call info8(1,0,0,
     .    ' sum ec : %;14,4D (SR) %;14,4D (FR) diff %;14,4D%N'//
     .    ' sum tc : %;14,4D (SR) %;14,4D (FR) diff %;14,4D',
     .    xv(1),sumecD,sumecD-xv(1),
     .    xv(2),sumtcD,sumtcD-xv(2),7,8)
      endif

      if (isw2 /= 0) then
C       call prrmsh('rho(SR)',rofi,rho,nr,nr,nsp)
        call dpzero(rho,nr*nsp)
        call dcopy(nr,rhoD,1,rho,1)
        call daxpy(nr,1d0,rhoD(1,2),1,rho(1,nsp),1)
C       call prrmsh('rho(D)',rofi,rho,nr,nr,nsp)
C       call rx('relativistic core still under development')
        sumec = 0 ; sumec(1) = sumecD
        sumtc = 0 ; sumtc(1) = sumtcD
      endif

      if (.not. associated(vloc,target=v)) deallocate(vloc)
      if (isw2 /= 0) deallocate(rhoD)

      end

      subroutine xyrhsr(ecore,l,z,nr,nre,g,rofi,wt,v,rho,deg,vrho,rhormx)
C- Make scalar Dirac density and integrate potential*density for one core state
C ----------------------------------------------------------------------
Ci Inputs
Ci   ecore :core eigenvalue
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   nr    :number of radial mesh points
Ci   nre   :Make density to the smaller of nr and nre
Ci   g     :normalized wave function times r
Ci   rofi  :radial mesh points
Ci   wt    :weights for numerical integration on radial mesh
Ci   v     :electronic part of spherical potential
Ci   deg   :occupation number
Co Outputs
Co   rho   :contribution core density from this state from 1..min(nr,nre)
Co   vrho  :integral
Co   rhormx:contribution to true density at nr from this state
Cr Remarks
Cr   xyrhsr makes the density at points 1..min(nr,nre), and the integral
Cr   of rho*density from 1..nre.  Thus:
Cr     if nre == nr, the density and integral are the same
Cr     if nre < nr  g**2 is negligible at nr (deep core state)
Cr     if nre > nr  (large sphere) contribution to rho*v is included
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh quadrature (radwgt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nre,l
      double precision z,ecore,vrho,rhormx
      double precision g(nr,2),rofi(nr),wt(nr),v(nr),rho(nr)
C ... Local parameters
      integer nrmx,ir
      double precision fpi,fllp1,r,tmc,c,gfac,rhoir,deg,rmax
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      fpi = 16d0*datan(1d0)
      fllp1 = l*(l+1)
      vrho = 0
      nrmx = min0(nr,nre)
C ... Make rho, and integrate vrho for points 1..nrmx
      do  ir = 2, nrmx
        r = rofi(ir)
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
        rho(ir) = rho(ir) + deg*rhoir
      enddo
      rmax = rofi(nr)
      rhormx = 0
C     nfp has the following line:
C     if (nre >= nr) rhormx = deg*g(nr,1)**2/(rmax*rmax*fpi)
      if (nre >= nr) rhormx = deg*rhoir/(fpi*rmax**2)

C ... Integrate rho*v from nrmx+1 .. nre
      do  ir = nrmx+1, nre
        r = rofi(ir)
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
      enddo

C     print *, 'xyrhsr : ', ddot(nre,rho,1,wt,1), vrho

      end
