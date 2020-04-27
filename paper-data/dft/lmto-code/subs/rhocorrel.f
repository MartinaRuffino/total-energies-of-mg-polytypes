      subroutine rhocorrel(isw,z,lmax,nsp,konfig,a,b,nr,rofi,v,g,
     .            kcor,lcor,qcor,tol,ec,sumec,sumtc,rho,gcore,ipr,
     .            ksop,nl,avw)
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
Ci   gn,fn  :work array holding normalized wave function times r (large and small components)
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   tol   :wave function tolerance
Ci   ipr   :0 no printout
Ci         :1 summary printout
Ci         :2 detailed printout
Ci   ksop  : spin-orbit coupling parameters
Ci   nl    :(global maximum l) + 1
Ci   avw   : average WS radius
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
Cr Remarks:
Cu Updates
Cu   28 Jan 13 core can be calculated with spin-averaged potential
Cu   19 Apr 01 core wave functions may be saved in gcore
C  ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer isw,ipr,lmax,konfig(0:lmax),nr,nsp,kcor,lcor,kcc,nl
      double precision a,b,z,g(nr,2),rho(nr,nsp),rofi(nr),tol,qcor(2),
     .  sumec(nsp),sumtc(nsp),gcore(nr,2,*),ec(*),avw
      real(8), target  :: v(nr,nsp)
      double precision ksop(0:nl-1,nsp,nsp,6)
      double precision pprel(5,0:lmax,2*(lmax+1),2,2)
C ... Local parameters
      character*1 pqn(9),ang(6)
      integer icore,isp,konf,l,ll,nodes,nre,stdo,intopt,isw0,isw1,imu,i
      double precision deg,dlml,e1,e2,ecor0,ecore,qcore,rhorim,rmax,
     .  rorim,slo,sum,tcore,val,vrho,kappa2,wt(nr),scalede
      real(8), pointer  :: vloc(:,:)
      double precision gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),
     . gsmt(2,2), Eout((lmax+1),(2*lmax+2))
      double precision gn(2,2,nr), fn(2,2,nr), gndot(2,2,nr),
     . fndot(2,2,nr),g2dot(2,2,nr), f2dot(2,2,nr)
      procedure(integer) :: nglob
C ... Heap
      data pqn /'1','2','3','4','5','6','7','8','9'/
      data ang /'s','p','d','f','g','h'/

      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      isw0 = mod(isw,10)
      isw1 = mod(isw/10,10)

      if (isw1 == 0) then
        vloc => v
      else
        allocate(vloc(nr,nsp))
        call dpscop(v,vloc,nr,1,1,0.5d0)
        call daxpy(nr,0.5d0,v(1,nsp),1,vloc,1)
        if (nsp == 2) then
          call dcopy(nr,vloc,1,vloc(1,2),1)
        endif
      endif

      rmax  = rofi(nr)
      b     = rmax/(dexp(a*nr-a)-1d0)
      call radwgt(intopt,rofi(nr),a,nr,wt)
      e1    = -2.5d0*z*z - 5d0
      e2    = 20d0
      icore = 0
!      do  isp = 1, nsp
      do  isp = 1, 1
        sumec(isp) = 0d0
        sumtc(isp) = 0d0
        qcore = 0d0
        rhorim = 0d0
        if (ipr >= 2) write(stdo,1)
        if (ipr >= 2.and.nsp == 2) write(stdo,'('' spin'',i2,'':'')')isp
    1   format(/' state  chg          ecor0',10x,'ecore',10x,'tcore',
     .     4x,'nre',2x,'rho(rmax)')
        do  l = 0, lmax
          do imu = 1,(2*l+2)
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
     .         g,sum,a,b,rofi,nr,nre,kcc)
             scalede = abs(ecor0)*1e-3
             print*,'E after rseq in rhocorrel =',ecor0!,'g(nr)=',
!     .  g(nr,1),'f(nr) =',g(nr,2)
             call rdeqcore(0,2*ecor0,ksop,z,v,ecor0,kcc,rofi,nr,nsp,a,
     .         b,l,(nl-1),nodes,imu,scalede,avw,g,gn,fn)
!             print*,'g11(nr)=',gn(1,1,nr),'f11(nr)=',fn(1,1,nr)
!             do i=1,nr
!                gn(:,:,i) = gn(:,:,i)*rofi(i)
!                fn(:,:,i) = fn(:,:,i)*rofi(i)
!             enddo
             ecore = ecor0
C     ... Correct core energy by using hankel bc's
             kappa2 = ecor0 - vloc(nr,isp) + 2*z/rmax
             if (nre == nr .and. kappa2 < 0d0) then
               dlml = -1d0 - dsqrt(-kappa2)*rmax
                 do  ll = 1, l
                   dlml = -kappa2*rmax*rmax/dlml - (2*ll+1)
                 enddo
               slo = val*(dlml+l+1)/rmax
               call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,vloc(1,isp),
     .           g,sum,a,b,rofi,nr,nre,kcc)
               print*,'RHOCORREL erseq = ',ecore
               call rdeqcore(0,2*ecore,ksop,z,v,ecore,kcc,rofi,nr,nsp,
     .            a,b,l,(nl-1),nodes,imu,scalede,avw,g,gn,fn)
!                  do i=1,nr
!                     gn(:,:,i) = gn(:,:,i)*rofi(i)
!                     fn(:,:,i) = fn(:,:,i)*rofi(i)
!                  enddo
             endif
             ec(icore) = ecore
             if (isw0 == 1) call dcopy(2*nr,g,1,gcore(1,1,icore),1)

C     --- Add to rho, make integral v*rho ---
!             do i=1,nr
!                    gn(:,:,i) = gn(:,:,i)*rofi(i)
!               fn(:,:,i) = fn(:,:,i)*rofi(i)
!             enddo
             call xyrhsrrel(ecore,l,z,nr,nre,gn,fn,rofi,wt,
     .         vloc(1,isp),rho(1,isp),deg,vrho,rorim)
             rhorim = rhorim + rorim
             tcore = ecore - vrho
             qcore = qcore + deg
             sumec(isp) = sumec(isp) + deg*ecore
             sumtc(isp) = sumtc(isp) + deg*tcore
               if (ipr >= 2) write (stdo,2) pqn(konf+1),ang(l+1),deg,
     .                                ecor0,ecore,tcore,nre,rorim
    2        format(1x,2a1,f8.2,3f15.6,i7,f9.5)
             enddo
          enddo
        enddo
        if (ipr > 0) write(stdo,3) qcore,sumec(isp),sumtc(isp),rhorim
    3   format(' sum q=',f5.2,'  sum ec=',f15.8,'  sum tc=',f15.8,
     .    '  rho(rmax)',f8.5)
      enddo
      end

      subroutine xyrhsrrel(ecore,l,z,nr,nre,gn,fn,rofi,wt,v,rho,deg,
     .  vrho,rhormx)
C- Make density and integrate potential*density for one core state
C ----------------------------------------------------------------------
Ci Inputs
Ci   ecore :core eigenvalue
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   nr    :number of radial mesh points
Ci   nre   :Make density to the smaller of nr and nre
Ci   gn,fn :normalized wave function times r (large and small components)
Ci   rofi  :radial mesh points
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
      double precision z,ecore,rhormx
      double precision gn(2,2,nr),fn(2,2,nr),rofi(nr),wt(nr),v(nr),
     .                 rho(nr)
C ... Local parameters
      integer nrmx,ir
      double precision fpi,fllp1,vrho,r,tmc,c,gfac,rhoir,deg,rmax
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
        rhoir = gfac*(gn(1,1,ir)**2+gn(2,2,ir)**2)
     .         + fn(1,1,ir)**2+fn(2,2,ir)**2
        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
        rho(ir) = rho(ir) + deg*rhoir/(2*l+2)
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
!        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        rhoir = gfac*(gn(1,1,ir)**2+gn(2,2,ir)**2)
     .         + fn(1,1,ir)**2+fn(2,2,ir)**2
        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
      enddo
C     vrho = vrho

      end

