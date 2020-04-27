      subroutine rnatm(pl,ql,n0,irchan,lmax,z,a,b,rofi,ev,
     .                 nr,rc,nsp,v,rho)
C- Renormalise charge or potential for a free atom
C ----------------------------------------------------------------------
Ci Inputs
Ci   pl    :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pl = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   ql    :valence charges for each l channel
Ci   n0    :dimensioning parameter
Ci   irchan:irchan(l+1) => suppress renormalization of that l
Ci   lmax  :maximum l for a given site
Ci   z     :nuclear charge
Ci   a     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rofi  :radial mesh points
Ci   ev    :valence eigenvalues
Ci   nr    :number of radial mesh points
Ci   rc    :Fermi radius for renormalization cutoff
Ci         :rc>0 => charge renormalization
Ci         :rc<0 => potential renormalization
Ci         :rc(2) is (optional) width.  If zero, choose rc/8.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cio Inputs/Outputs
Cio  v     :On input, spherical potential (atomsr.f)
Cio        :if rc<0, renormalized on output  NOT IMPLEMENTED
Co Outputs
Co   rho   :renormalized on output
Cl Local variables
Cl   ltop  :max l for which charge is nonzero
Cr Remarks
Cr
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   01 Feb 06 Adapted from mol/rnatm.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,n0,lmax,irchan(1+lmax)
      double precision z,a,b,rc(2),pl(n0,2),ql(n0,2)
      double precision rho(nr,nsp),v(nr,nsp),rofi(nr),ev(25)
C ... Local parameters
      logical pot
      integer intopt,iprint,ir,isp,konfig,l,lp1,ltop,stdo
      double precision sum1(10,2),sum2(10,2),tol,r,fac,sum,rmax,rcw
      double precision rh1(nr),rh2(nr),rwgt(nr),ddot
      procedure(integer) nglob
C ... External calls
      external dpcopy,info2,makrvl,radwgt,rseq,rx

C     call prrmsh('input rho',rofi,rho,nr,nr,nsp)

      if (rc(1) == 0) return

C --- Setup ---
      stdo = nglob('stdo')
      rmax = rofi(nr)
      intopt = 10*nglob('lrquad')
      call radwgt(intopt,rmax,a,nr,rwgt)
      rcw = rc(2)
      if (rcw == 0) rcw = abs(rc(1)/8)
      tol = 1d-10
      if (lmax > 8) call rx('RNATM: lmax > 8')
      call info2(20,1,0,
     .  ' RNATM: renormalize sphere density  rc=%,1d  w=%,1d',rc,rcw)

      pot = rc(1) < 0
      if (pot) call rx('RNATM not set up for pot now')

      ltop = -1
      do   isp = 1, nsp
        if (nsp == 2) call info2(20,0,0,' Spin %i:',isp,0)
        do  l = 0, lmax
        lp1 = l+1
        konfig = pl(lp1,isp)
        sum1(lp1,isp) = 0
        sum2(lp1,isp) = 0
        if (ql(lp1,isp) < 1d-6) goto 20
        call makrvl(z,l,a,b,nr,rofi,konfig,v(1,isp),ev,tol,rh1)
C       call intrho(rh1,rofi,a,b,nr,sum1(lp1,isp))
        sum1(lp1,isp) = ql(lp1,isp) * ddot(nr,rwgt,1,rh1,1)
        if (irchan(lp1) == 0) then
          do  ir = 1, nr
            r = rofi(ir)
            fac = 1d0/(dexp((r-dabs(rc(1)))/rcw) + 1)
            rh2(ir) = rh1(ir)*fac
          enddo
        else
          call dpcopy(rh1,rh2,1,nr,1d0)
        endif
C       call intrho(rh2,rofi,a,b,nr,sum2(lp1,isp))
        sum2(lp1,isp) = ql(lp1,isp) * ddot(nr,rwgt,1,rh2,1)
        if (ql(lp1,isp) < 1d-6) goto 20
        fac = sum1(lp1,isp)/sum2(lp1,isp)
        do  ir = 1, nr
          rho(ir,isp) = rho(ir,isp) + ql(lp1,isp)*(fac*rh2(ir)-rh1(ir))
        enddo
        ltop = max(ltop,l)
   20   continue
        enddo

C --- Printout ---
      if (iprint() >= 20) then
        sum = ddot(nr,rwgt,1,rho(1,isp),1)
        write(stdo,333) (sum1(lp1,isp), lp1=1,ltop+1)
        write(stdo,334) (sum2(lp1,isp), lp1=1,ltop+1)
        write(stdo,335) sum
  333   format(' Ql before scaling:',8f9.5)
  334   format(' Ql after  scaling:',8f9.5)
  335   format(' valence q after renormalisation:',f9.5)
      endif
      enddo

      end
      subroutine makrvl(z,l,a,b,nr,rofi,konfig,v,ev,tol,rho)
C- Makes contribution to charge density from a given spin and l,
C  for free atom and given spherical potential (Adapted from nwrofp)
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :complex energy
Ci   z     :nuclear charge
Ci   l     :charge for particular l
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   konfig:core configuration
Ci   v     :spherical potential (atomsr.f)
Ci   ev
Ci   tol
Co Outputs
Co   rho  :partial density for given l
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   01 Feb 06
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ir,ival,konfig,l,nn,nr,nre,jr
      double precision v(nr),rofi(nr),rho(*),ev(*)
      double precision a,b,z,tol
C ... Local parameters
      integer lrel,kc,nrx
      parameter( nrx=5001 )
      double precision eb1,eb2,eval,slo,sum,val,gfac,fllp1,c,tmc,r
      double precision g(2*nrx)
      procedure(integer) nglob
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

C     stdo = nglob('stdo')
      eb1 = -50.d0
      eb2 = 15d0
      ival = 0
      nn = konfig-(l+1)
      ival = ival+1
      eval = ev(ival)
      val = 1.d-30
      slo = -val
      lrel = mod(nglob('lrel'),10)
      if (lrel == 0) then
        call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v,g,
     .    sum,a,b,rofi,nr,nre,kc)
        do  77  ir = 1, nre
   77   rho(ir)= a*(rofi(ir)+b)*g(ir)**2
      else
        call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v,g,
     .    sum,a,b,rofi,nr,nre,kc)
        fllp1 = l*(l+1)
c       c = 274.074d0
        rho(1) = 0
        do  78  ir = 2, nre
          jr = ir+nr
          r = rofi(ir)
          tmc = c - (v(ir) - 2*z/r - eval)/c
          gfac = 1 + fllp1/(tmc*r)**2
          rho(ir) = gfac*g(ir)**2 + g(jr)**2
   78   continue
      endif
      do  79  ir = nre+1, nr
   79 rho(ir)= 0d0

C      print *, 'makrvl: l=',l
C      call prrmsh('psi**2 in makrvl',rofi,rho,nre,nre,1)

      end
