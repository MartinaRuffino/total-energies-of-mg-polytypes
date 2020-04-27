      subroutine rdeqmu(enu,gfn,iopt,z,v,rofi,nr,nsp,kc,a,b,l,lmx,imu,psi,psd,gsmt,GFr)
C- Solves radial Dirac equations and energy derivatives for a given potential, l, mu
C ----------------------------------------------------------------------
Ci Inputs
Ci   enu   :linearization energy.
Ci         :enu is normally an input, but if iopt>1, it can be output
Ci   iopt  :1s digit
Ci         :0 linearization energy (input).  The full Dirac 4-vector is calculated.
Ci         :1 input E is an estimate only; its value is determined internally
Ci         :  as a function of the logarithmic derivative of the dominant component,
Ci         :  specified through gsmt, only the log derivative of the dominant
Ci         :  large component is matched, providing an incompletely determined
Ci         :  solution of the Dirac equation.
Ci         :  You must supply gsmt to use this mode.
Ci         :2 similar to iopt=1, except that the cross-coupling between (g1,f1) and
Ci         :  (g2,f2) through the magnetic field is omitted.  This ensures that
Ci         :  (g2,f2) vanish for alpha=1 and (g1,f1) vanish for alpha=2; the result
Ci         :  satisfies an approximate differential equation
Ci         :  You must supply gfn to use this mode.
Ci         :10s digit
Ci         :  0 initialize gfn to 0; compatible only with 1s digit iopt = 0
Ci         :  1 use given gsmt to populate gfn
Ci         :    For 1s digit=0, this is an initial estimate; it is updated
Ci         :    For 1s digit>0, it supplies a fixed boundary condition.
Ci         :100s digit
Ci         :  0 return psi, psd, GFr
Ci         :  1 return psi only
Ci   z     :nuclear charge
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   kc    :rofi(kc) is point where inwards and outwards solutions are joined
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l     :l quantum number
Ci   lmx   :dimensions GFr.  lmx<l => do not return GFr
Ci   imu   :projection of j (positive integer, counts -l-1/2...l+1/2)
Cio Inputs/Outputs
Cio  gfn   :Boundary conditions for Dirac equation, 8 components to include alpha=1,2; see rdeq0
Cio        :gfn(1:2,:,:,:) large/small component
Cio        :gfn(:,1:2,:,:) g1, g2
Cio        :gfn(:,:,1:2,:) for psi or psid (psid is used only to improve convergence)
Cio        :gfn(:,:,:,1:2) for alpha = 1,2
Co         :Input gfn: usage depends on 10s digit iopt
Co         : 0 gfn is initialized to zero (rdeq finds values internally)
Co         :>0 Use as initial boundary condition for rdeq0
Co         :Output gfn:
Co         :returned by rdeq0.  How gfn is modified depends on choice of b.c. (see 1s digit iopt)
Co Outputs
Co   psi   :psi(meshpoint,1,N,alpha) large component r*g, called P in Ebert's paper
Co         :psi(meshpoint,2,N,alpha) small component r*f, called Q in Ebert's paper
Co         :N = 1 or 2 corresponding to kappa (kappa = l or -l-1)
Co         :alpha = 1 or 2, called lambda in Schick's paper (PhysRevB, 54 1610 (1996))
Co         :                called i in Ebert's paper (J Phys Condens. Matter 1, 9111 (1989)
Co         :psi(:,1:2,1:2,alpha) must be solved as four coupled equations
Co   psd   :energy derivatives of psi
Co   gsmt  :gsmt(N,alpha) = slope dg/dr of g_N,alpha at augmentation radius, N=1,2 alpha=1,2
Co   GFr   :GFr(4,4,2,2,l,imu,4):decomposition of overlap integrals for the calculation of averages [Turek (6.141)]
Co           1) <g|g>,<f|f>, 2) <g|gdot>,<f|fdot>,3) <gdot|g> ... 4) <gdot|gdot>,...
Cl Local variables
Cr Remarks
Cr   REWRITE ... equations are wrong ... see rdeqinout
Cr   The large component g and small component f in the Dirac equation
Cr   for a central potential V(r) obey the following equation (B=0)
Cr       (d/dr + (1+k)/r) g_k(E,r)  - [1+(E-V)/c^2]cf_k(E,r) = 0     (N1)
Cr       (d/dr + (1-k)/r) cf_k(E,r) + [(E-V)/c^2]   g_k(E,r) = 0     (N2)
Cr   k refers to kappa: kappa = l corresponding to j=l-1/2,  or to -l-1 (j=l+1/2)
Cr   and c to the speed of light.
Cr   r*g and r*f are regular for r -> 0 (both vanish at r=0).
Cr   Note that these equations do not depend on the quantum number mu since B=0.
Cr
Cr   (g,f) are normalized within the Wigner-Seitz sphere such that
Cr     int[g**2(E,r) + f**2(E,r)] r^2 dr = 1
Cr
Cr   For small r, the potential is dominated by the nuclear 2*Z/r.
Cr   Choose a trial solution  g(r) := r^a ;  f(r) := r^a * q;
Cr   Eq. (N2) yields a relation
Cr     q = (c/2Z) (kappa + a+1) = (c/2Z) (kappa + gamma)  defining gamma=a+1
Cr
Cr   and the solution can be expressed as :
Cr       g(r) = r^(gamma-1) Sum_n=1^infty (p_n r^n)
Cr       f(r) = r^(gamma-1) Sum_n=1^infty (q_n r^n)
Cr   where
Cr     gamma = sqrt[kappa^2 - (2*Z/c)^2]
Cr     q/p   = (c/2Z) (kappa + gamma)
Cr   When a magnetic field is present there is an additional term
Cr      xxx
Cr   In that case there are two independent solutions labelled by alpha=1,2
Cr   with different boundary conditions.  For r->0, these conditions are:
Cr     alpha=1 => gk1 and fk1 are proportional to r^gam for small r
Cr                gk2 and fk2 are zero
Cr     alpha=2 => gk2 and fk2 are proportional to r^gam for small r
Cr                gk1 and fk1 are zero
Cr   In the absence of a magnetic field the two solutions are equivalent with
Cr   the transposition of alpha=1,2.  Only two of the components of each solution
Cr   are nonzero.
Cr
Cr   In this version enu depends only on kappa
Cu Updates
Cu   01 Aug 16 (MvS) Redesign of solver
Cu   30 Jul 16 (MvS) First attempt at redesign of solver
Cu   24 Oct 14 (Vishina) Rewrite of rseq.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,imu,l,lmx,kc,iopt
      double precision enu,z,a,b
      double precision gfn(2,2,2,2),v(nr,nsp),rofi(nr,2),gsmt(2,2,2)
      double precision psi(nr,2,2,2),psd(nr,2,2,2)
C ... Local parameters
      integer i,alp,a1,a2,k1,k2,idmod
      double precision gfx(2,2,2,2)
      double precision psie(nr,2,2,2,4)
      double precision ei(4),g0(2),g0x(2)
      double precision dele,de1,de2,del1,del2,errmx
      double precision gsmt1(2,2)
      double precision GFr(4,4,2,2,0:lmx,2*(lmx+1),4)
      double precision vm(nr+2),bm(nr+2) ! Spin averaged potential and Magnetic field
      procedure(real(8)) :: dot3,diracp

C     idmod>0 => freeze f/g rmt; idmod=0 => freeze enu
      idmod = mod(iopt,10)

C     Split potential into spin-averaged part + magnetic field;
C     Extrapolate by 2 radial mesh points to allow Runge Kutta integration to nr
      vm(1:nr) = (v(1:nr,1) + v(1:nr,nsp))/2
      call vpxtrap(nr,2,5,a,b,rofi,vm,errmx)
      if (nsp == 1) then
        bm(1:nr) = 0
      else
        bm(1:nr) = ((v(1:nr,2) - v(1:nr,1))/2)
      endif
      call vpxtrap(nr,2,5,a,b,rofi,bm,errmx)

C --- Make normalized scalar g, and energy derivatives at MT boundary
C      vm(:) = (v(:,1) + v(:,2))/2
C      dele = .002d0; ii(1:3) = [-1, 1, 0]
C      do  j1 = 1, 3
C        j2 = ii(j1)
C        call rsq1(0,enu+j2*dele,l,z,vm,nr,g,val,slo,i,a,b,rofi,nr)
C        call prodrw(nr,2,g,g,rofi(1,2),norm(j2)); norm(j2) = sqrt(norm(j2))
C        call dscal(2*nr,1/norm(j2),g,1)
C        if (j1 == 1) then
C          gdotn = g(nr)
C          fdotn = g(2*nr)
C        elseif (j1 == 2) then
C          gdotn = (g(nr) - gdotn)/(2*dele)
C          fdotn = (g(2*nr) - fdotn)/(2*dele)
C        endif
C      enddo
CC     norm = norm/product
CC     Starting boundary conditions
C      g0 =  0                       ! Initial condition unknown
C      gfn = 0
C      gfn(1,1,1,1) = g(nr)            ! large component from scalar rel, alpha=1
C      gfn(2,1,1,1) = g(2*nr)          ! small component from scalar rel, alpha=1
C      gfn(1,2,1,2) = g(nr)            ! large component from scalar rel, alpha=2
C      gfn(2,2,1,2) = g(2*nr)          ! small component from scalar rel, alpha=2
C
CC      gfn(1,1,2,1) = ()/(2*dele)
CC      gfn(2,1,2,1) = (slo(1)/norm(1)-slo(-1)/norm(-1))/(2*dele)
CC      gfn(1,2,2,2) = (val(1)/norm(1)-val(-1)/norm(-1))/(2*dele)
CC      gfn(2,2,2,2) = (slo(1)/norm(1)-slo(-1)/norm(-1))/(2*dele)
C
C      gfn(1,1,2,1) = gdotn
C      gfn(2,1,2,1) = fdotn
C      gfn(1,2,2,2) = gdotn
C      gfn(2,2,2,2) = fdotn

C --- Dirac equation ---
      g0 =  0                       ! r->0 boundary condition unknown
      if (mod(iopt,10) == 0) then   ! Otherwise, gfn is given
        gfn = 0
      endif
      call rdeq0(enu,z,vm,bm,rofi,nr,kc,a,b,l,imu,idmod,g0,gfn,psi,gsmt)

      if (iopt >= 100) return

C --- Energy derivative of psi (5 point numerical differentiation) ---
C     print *, '!! skip e derivatives'; gdot = 0; fdot = 0; return
      dele = .002d0
      ei(1) = 1
      ei(2) = -1
      ei(3) = 1.5d0
      ei(4) = -1.5d0
      do  i = 1, 4
        g0x = g0
        ei(i) = enu + dele*ei(i)
        gfx = gfn
C       this is a good guess, but seems to be more efficient to internally remake gfn
C        do  j1 = 1, 2
C          gfx(:,j1,1,j1) = gfn(:,j1,1,j1) - dele*ei(i) * gfn(:,j1,2,j1)
C        enddo
        call rdeq0(ei(i),z,vm,bm,rofi,nr,kc,a,b,l,imu,0,g0x,gfx,psie(1,1,1,1,i),gsmt1)
      enddo
      de1  = (ei(1) - ei(2))/2
      del1 = (ei(1) + ei(2))/2 - enu
      de2  = (ei(3) - ei(4))/2
      del2 = (ei(3) + ei(4))/2 - enu
      do  a1 = 1, 2
        do  a2 = 1, 2
          call dfphi(de1,del1,de2,del2,nr,psi(1,1,a1,a2),psie(1:nr,1,a1,a2,1:4),.true.)
          call dfphi(de1,del1,de2,del2,nr,psi(1,2,a1,a2),psie(1:nr,2,a1,a2,1:4),.true.)
          psd(1:nr,1,a1,a2) = psie(1:nr,1,a1,a2,1)
          psd(1:nr,2,a1,a2) = psie(1:nr,2,a1,a2,1)
        enddo
      enddo
C      call prrmsh('g11-dot',rofi,psd(1:nr,1,1,1),nr,nr,1)
C      call prrmsh('f22-dot',rofi,psd(1:nr,2,2,2),nr,nr,1)

C     Debugging
C      call info5(70/10,0,0,' <psii psij>i l=%i  imu=%i  %3;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psi,psi,rofi(1,2)),
C     .   diracp(nr,psi,psi(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psi(1,1,1,2),rofi(1,2))],4,5)
C      call info5(70/10,0,0,' <psii dotj>i l=%i  imu=%i  %4;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psi,psd,rofi(1,2)),
C     .   diracp(nr,psi,psd(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psd,rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psd(1,1,1,2),rofi(1,2))],4,5)
C      call info5(70/10,0,0,' <doti dotj>i l=%i  imu=%i  %3;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psd,psd,rofi(1,2)),
C     .   diracp(nr,psd,psd(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psd(1,1,1,2),psd(1,1,1,2),rofi(1,2))],4,5)

C ... Render psd(alpha=2) and psd(alpha=1) orthogonal
C     Don't use: self-consistent moment is too small in Fe
C     However, if allowed to be nonorthogonal charge neutrality is not quite conserved.
C     Remedy: set qnur(2,:,:,1,2) and qnur(2,:,:,2,1) to zero (see qrel2z12)
C      if (imu /= 1 .and. imu /= 2*l+2) then
CC       Orthogonalize psd(alpha=2) to psd(alpha=1)
C        call orthdirac(1,nr,rofi(1,2),psd,psd(1,1,1,2))
CC       Symmetrically orthogonalize psd(alpha=1) and psd(alpha=2)
CC       call orthdirac(2,nr,rofi(1,2),psd,psd(1,1,1,2))
C      endif

C ... Orthogonalize phidot to phi
      do  alp = 1, 2
        call orthdirac(1,nr,rofi(1,2),psi,psd(1,1,1,alp))
        call orthdirac(1,nr,rofi(1,2),psi(1,1,1,2),psd(1,1,1,alp))
      enddo

C     Debugging
C      call info5(70/10,0,0,' <psii psij>f l=%i  imu=%i  %3;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psi,psi,rofi(1,2)),
C     .   diracp(nr,psi,psi(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psi(1,1,1,2),rofi(1,2))],4,5)
C      call info5(70/10,0,0,' <psii dotj>f l=%i  imu=%i  %4;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psi,psd,rofi(1,2)),
C     .   diracp(nr,psi,psd(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psd,rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psd(1,1,1,2),rofi(1,2))],4,5)
C      call info5(70/10,0,0,' <doti dotj>f l=%i  imu=%i  %3;14,9D',
C     .  l,imu,
C     .  [diracp(nr,psd,psd,rofi(1,2)),
C     .   diracp(nr,psd,psd(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psd(1,1,1,2),psd(1,1,1,2),rofi(1,2))],4,5)

      if (l > lmx) return

C ... Resolve overlap <psi_i | psi_j> into the four components
      do  a1 = 1, 2
        do  a2 = 1, 2
          do  k1 = 1, 4
            do  k2 = 1, 4
              if (mod(k1,2) == mod(k2,2)) then ! don't need <g|f>
                GFr(k1,k2,a1,a2,l,imu,1) = dot3(nr,psi(1,k1,1,a1),psi(1,k2,1,a2),rofi(1,2))
                GFr(k1,k2,a1,a2,l,imu,2) = dot3(nr,psi(1,k1,1,a1),psd(1,k2,1,a2),rofi(1,2))
                GFr(k1,k2,a1,a2,l,imu,3) = dot3(nr,psd(1,k1,1,a1),psi(1,k2,1,a2),rofi(1,2))
                GFr(k1,k2,a1,a2,l,imu,4) = dot3(nr,psd(1,k1,1,a1),psd(1,k2,1,a2),rofi(1,2))

C                if (a1 == a2)
C     .            call info5(10,0,0,' GFr a = %2,2i  k = %2,2i  l imu = %2,2i %4;12,7D',
C     .            [a1,a2],[k1,k2],[l,imu],GFr(k1,k2,a1,a2,l,imu,1:4),5)
              endif
            enddo
          enddo
        enddo
      enddo

C     call prmx('GFr [Turek (6.141)]',GFr,16,16,size(GFr)/16)

      end subroutine rdeqmu

      subroutine rdeqcore(iopt,E,z,v,rofi,nr,nsp,nn,kc,a,b,l,Ekmu,rhoD,tcorD,nre,rhormxD)
C- Core eigenvalues and charge density from the Dirac equation, for a given l
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :1 Full Dirac equation
Ci         :2 Dirac equation neglecting coupling between lambda=1 and lambda=2
Ci         :  in which case 4 vector simplies to 2-vector
Ci         :10s digit
Ci         :0 Do not accumulate the density rhoD
Ci         :1 Add to rhoD for each of the 2l+2 core levels
Ci         :2 Same as 1, but rhoD is initialized on input
Ci   E     :trial eigenvalue (from scalar Dirac)
Ci   z     :nuclear charge
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nn    :number of nodes (used to determine sign of psi(r->0)
Ci   kc    :Turning point where inward and outward solutions are joined
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l     :orbital angular momentum
Co Outputs
Co   Ekmu  :Eigenvalues of Dirac equation
Co   psi   :wave function (not returned now)
Co   rhoD  :core density times 4*pi*r*r
Cl Local variables
Cl         :
Cr Remarks
Cr  Approximately adapt Ebert's method (J. Phys. Cond. Matt. 1, 9111 (1989)
Cr  psi consists of 4 components (g,f)_lam, lam = 1, 2.
Cr  psi is built out of a combination of two independent solutions alpha=1,2
Cr  psi(1,lam,alpha) large component g, lam=1,2
Cr  psi(2,lam,alpha) small component f, lam=1,2
Cr  Each of the alpha=1,2 solns can be integrated to rofi(kc) from rmax inwards,
Cr  or from the origin outwards
Cr  Three independent parameters are available to assemble a solution:
Cr    ratio of outgoing alpha=2/alpha=1 solutions
Cr    ratio of incoming alpha=2/alpha=1 solutions
Cr    ratio of incoming to outgoing solution
Cr  Energy supplies fourth degree of freedom to satisfy the four matching conditions.
Cr
Cr  Let gao = g1 component of outward solution alpha (dominant for alpha=1: kap=l unless imu=1 or 2l+2)
Cr  Let fao = f1 component of outward solution alpha
Cr  Let Gao = g2 component of outward solution alpha (dominant for alpha=2: kap=-l-1)
Cr  Let Fao = f2 component of outward solution alpha
Cr
Cr  In either case write solution as follows : r<rofi(kc)
Cr        r < rofi(kc)       r > rofi(kc)
Cr    go = g1o + Xo g2o   gi = A g1i + Xi g2i
Cr    fo = f1o + Xo f2o   fi = A f1i + Xi f2i
Cr    Go = G1o + Xo G2o   Gi = A G1i + Xi G2i
Cr    Fo = F1o + Xo F2o   Fi = A F1i + Xi F2i
Cr
Cr  Matching conditions : determine Xo, Xi, A, E that satisfy:
Cr    g1o + Xo g2o - A g1i - Xi g2i = 0
Cr    f1o + Xo f2o - A f1i - Xi f2i = 0
Cr    G1o + Xo G2o - A G1i - Xi G2i = 0
Cr    F1o + Xo F2o - A F1i - Xi F2i = 0
Cr
Cr  For a given energy this results in the following underdetermined linear equations
Cr    (g2o - g1i - g2i) (Xo) = -g1o    (M1)
Cr    (f2o - f1i - f2i) (A)  = -f1o    (M2)
Cr    (G2o - G1i - G2i) (Xi) = -G1o    (M3)
Cr    (F2o - F1i - F2i)      = -F1o    (M4)
Cr  Energy supplies a fourth degree of freedom to satisfy all four equations.
Cr
Cr  Practical method of solution:
Cr  Use the fact that M1 is dominant for lam=1 and M3 dominant for lam=2.
Cr  Solve as follows:
Cr  lam=1: M1,M3,M4 form a set of linear conditions on Xo, A, Xi, with A near 1.
Cr    Determine Xo, A, Xi from (M1,M3,M4)
Cr    Discontinuity in M2 approximately fixes discontinuity in slope of g1.
Cr    Estimate correction to energy and iterate to solution (similar to method in rdeq0)
Cr  lam=2: M1,M2,M3 form a set of linear conditions on Xo, A, Xi, with A near 1.
Cr    Determine Xo, A, Xi from (M1,M2,M3)
Cr    Discontinuity in M4 approximately fixes discontinuity in slope of g2.
Cr    Estimate correction to energy and iterate to solution (similar to method in rdeq0)
Cr
Cr  p core levels: In the nonmagnetic case, the 2x3 p orbitals split into :
Cr  From rdeqmu for small r : g(r) = r^(gamma-1) Sum_n=1^infty (p_n r^n)
Cr  where gamma = sqrt[kappa^2 - (2*Z/c)^2]
Cr  In the limit c->infty, g(r) ~ r^|kap| for small r
Cr     l  imu  kap  mu  lam  S
Cr     0   1   -1  -0.5      3/2
Cr     0   1   -1   0.5      3/2
Cr     l  imu  kap  mu  lam  P
Cr     1   1   -2  -1.5      3/2
Cr     1   2    1  -0.5  1   1/2 g ~ r
Cr     1   2   -2  -0.5  2   3/2
Cr     1   3    1   0.5  1   1/2 g ~ r
Cr     1   3   -2   0.5  2   3/2
Cr     1   4   -2   1.5      3/2
Cr     l  imu  kap  mu  lam  D
Cr     2   1   -3  -2.5      5/2
Cr     1   2    2  -1.5  1   3/2
Cr     2   2   -3  -1.5  2   5/2
Cr     2   3    2  -0.5  1   3/2
Cr     2   3   -3  -0.5  2   5/2
Cr     2   4    2   0.5  1   3/2
Cr     2   4   -3   0.5  2   5/2
Cr     2   5    2   1.5  1   3/2
Cr     2   5   -3   1.5  2   5/2
Cr     2   6   -3   2.5      5/2
Cr
Cb Bugs
Cb   calls info ... should not do this since may be called in parallel!
Cu Updates
Cu   15 Aug 16 First created ... makes wave functions only so far
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iopt,nr,nsp,l,kc,nn
      double precision z,a,b,E,tcorD,rhormxD,Ekmu(2*(2*l+1))
      double precision v(nr,nsp),rofi(nr,2),rhoD(nr,2)
C     double precision psi(nr,4,2*(2*l+1))
C ... Dynamically allocated arrays
      real(8), allocatable :: psi4o(:,:,:),psi4i(:,:,:),rofe(:,:),rhom(:),amom(:)
C ... Local parameters
      logical :: userfalsi,norfalsi
C     logical :: debug=.false.
      integer, parameter :: IPRT1=60, IPRT2=70, IPRT3=100, maxit=100
      integer i,kap,imu,ir,alp,lam,lamb,iter,ipr,k1,k2,jpr,nre,icore,nl,opt0,opt1
C     integer idx(2,(l+1),2*(l+1))
      double precision sre,mu,q,dq,c,E0,Ealp(2),errmx,enu,de,growth,wkg(28)
      double precision gsmt(2,2),sec(3,3),seci(3,3),rhs(3),Pv(2,2),g0(2),tcorm
      double precision gfn(2,2,2,2)      ! For Boundary conditions at rmax
      double precision vm(nr+2),bm(nr+2) ! Spin averaged potential and Magnetic field
      double precision valslo(2,2,2),vrhom,rhormxm
C     real(8) :: psi4o(nr,4,2),psi4i(nr,4,2)
      real(8), parameter :: tolq = 1d-6
      real(8), parameter :: tolgrmax=1d-25 ! wave function < tolgrmax is assumed 0
      real(8), parameter :: tolg=1d-13   ! Tolerance in w.f. discontinuity
      real(8), parameter :: tole=1d-13   ! Tolerance in energy

C     procedure(real(8)) :: dot3,diracp
      procedure(integer) :: isw,nglob
      procedure(real(8)) :: ddot,diracp

      data wkg /28*0d0/

C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      call getpr(ipr)
      call pshpr(ipr)  ! So the verbosity can be locally manipulated for debugging
      opt0 = mod(iopt,10)
      opt1 = mod(iopt/10,10)
      g0 = 0
      sre = dsqrt(dabs(-E-E**2/c/c))
      nl = l+1
      if (opt1 == 2) call dvset(rhoD,1,nr*2,0d0)
      tcorD = 0

C     print *, '!! relativistic core'; call setpr(70)

C     call mstokm(2,nl,1,nl*nl,sec,sec,idx)   ! Get idx for this l

C     Find nre beyond which w.f. can be assumed 0.  Roughly, psi = exp(-sre*r)
      nre = nr
      q = -dlog(tolgrmax)/sre
      if (q < rofi(nr,1)) then
        call huntx(rofi,nr,q,0,nre)
        nre = max0(51,((nre-1)/2)*2+1)
      endif
      allocate(psi4o(nre,4,2),psi4i(nre,4,2),rofe(nre,2),rhom(nre),amom(nre))
C     allocate(psi4o(nre,4,2),psi4i(nre,4,2*2),rofe(nre,2),rhom(nre),amom(nre))
      call dpzero(rhom,size(rhom))
      call dpzero(amom,size(rhom))
      call dcopy(nre,rofi,1,rofe,1)
      call radwgt(10*nglob('lrquad'),rofe(nre,1),a,nre,rofe(1,2))

      icore = 0 ! index to core levels.  Levels ordeed imu=1,2(lam=1),2(lam=2),...2*l+1
      gfn = 0
      gfn(1,1,1,1) = dexp(-sre*rofe(nre,1)) ! Estimate for psi(nre)
      if (mod(nn,2) /= 0) gfn(1,1,1,1) = -gfn(1,1,1,1) ! Odd number of nodes
      gfn(1,2,1,2) = gfn(1,1,1,1)

C     Split potential into spin-averaged part + magnetic field
C     Extrapolate by 2 radial mesh points to allow Runge Kutta integration to nre
      vm(1:nr) = (v(1:nr,1) + v(1:nr,nsp))/2
      if (nre > nr-2) call vpxtrap(nre,2,5,a,b,rofe,vm,errmx)
      if (nsp == 1) then
        bm(1:nr) = 0
      else
        bm(1:nr) = ((v(1:nr,2) - v(1:nr,1))/2)
      endif
      if (nre > nr-2) call vpxtrap(nre,2,5,a,b,rofe,bm,errmx)

C --- Accumulate density for each quantum state mu ---
      do  imu = 1, 2*l+2
        mu  = imu - l - 1.5d0
C       Vary boundary condition until it is consistent with E
        dq = 9999
        do  i = 1, 5
          sre = dsqrt(dabs(-E-E**2/c/c))
          q = -sre/(1d0+E/c/c)/c
          gfn(2,1,1,1) = q*gfn(1,1,1,1)
          gfn(2,2,1,2) = q*gfn(1,2,1,2)
          E0 = E; Ealp = E
C         if (imu == 2 .and. l == 1) call snot(1,debug)
          call rdeq0(Ealp,z,vm,bm,rofe,nre,kc,a,b,l,imu,opt0,g0,gfn,psi4i,gsmt)
C         call testdirac(nre,l,imu,Ealp(1),z,rofe,psi4i,vm)
          if (imu == 1 .or. imu == 2*l+2) Ealp(2) = Ealp(1)
          E = Ealp(1)
          dq = q + sre/(1d0+E/c/c)/c
          call info8(IPRT1,1,0,' rdeqc imu=%i E0=%;9F update q=%;9F : '//
     .      '%s,E=%?#(n)#%;11F#%2;11F# => dq=%1,3;3g%?#n==0# ... repeat# ... bc cnvgd#',
     .      imu,E0,q,isw(imu == 1.or.imu == 2*l+2),Ealp,dq,isw(abs(dq) < tolq),8)
          if (abs(dq) < tolq) exit
          if (i == 5) call rx('rdeqcore failed to find b.c. consistent with energy')
        enddo
C        call info5(1,0,0,' l=%i imu=%i psi1(kc) %4;18,12D %;18,12D',l,imu,psi4i(kc,:,1),
C     .    diracp(nre,psi4i(1,1,1),psi4i(1,1,1),rofe(1,2)),5)
C        call info5(1,0,0,' l=%i imu=%i psi2(kc) %4;18,12D %;18,12D',l,imu,psi4i(kc,:,2),
C     .    diracp(nre,psi4i(1,1,2),psi4i(1,1,2),rofe(1,2)),5)

C       call prrmsh('psi4 from rdeq0',rofe,psi4i,nre,nre,4)

C   ... Handle case that mu is extremal value => one solution
        if (imu == 1 .or. imu == 2*l+2) then

          call info5(IPRT1,0,0,' rdeqc cnvg l=%i, kap=%i, mu=%d %33pnod=%i  E=%;11F',l,-l-1,mu,nn,E)
          icore = icore+1
          Ekmu(icore) = E

C         Renormalize using quadrature containing mesh of all points 1:nr
          call dscal(4*nre,1/dsqrt(diracp(nre,psi4i,psi4i,rofi(1,2))),psi4i,1)
C         Make density and magnetization
          call xyrhfr(Ekmu(icore),l,imu,z,nre,nre,psi4i,rofe,rofe(1,2),vm,bm,rhom,amom,vrhom,tcorm,rhormxm)
C          if (debug) then
C            print *, 'imu=1 or 2*l+2, l=',l
C            call prrmsh('psi4i imu=1 or 2*l+2',rofe,psi4i,nre,nre,4)
C          endif

          tcorD = tcorD + tcorm
          rhormxD = rhormxD + rhormxm
          if (opt1 >= 1) then
            do  i = 1, 2
              call daxpy(nre,1d0/2,rhom,1,rhoD(1,i),1)
              call daxpy(nre,dble(2*i-3)/2,amom,1,rhoD(1,i),1)
            enddo
          endif

C       Find linear combination of (alpha=1, alpha=2) matching at kc
        else
C fplot -x 0,1 -ord 'x2' psi1  -ord 'x3*274' -lt 1,col=1,0,0 psi1   -lt 2,col=0,0,0 -ord 'x4*274' psi1   -lt 2,col=1,0,0 -ord 'x5*274^2' psi1
C fplot -x 0,1 -ord 'x4' psi2  -ord 'x5*274' -lt 1,col=1,0,0 psi2   -lt 2,col=0,0,0 -ord 'x2*274' psi2   -lt 2,col=1,0,0 -ord 'x3*274^2' psi2
C          call prrmsh('psi1',rofe,psi(1,1,1,1),nre,nre,4)
C          call prrmsh('psi2',rofe,psi(1,1,1,2),nre,nre,4)

C     ... Find solutions lam=1,2 corresponding to kap = l and kap = -l-1.
C         Only 1 solution if imu=1 or 2*l+2.
          do  lam = 1, 2
          kap = l; if (lam == 2) kap = -l-1
          enu = Ealp(lam)
          lamb = 3-lam

          ir = 0
          iter = 0; ir = 0; userfalsi = .false.; growth = 1
          norfalsi = .false.    ! When true, rfalsi acceleration is not used
          do  while (iter <= maxit)
          iter = iter+1

          do  alp = 1, 2
            call rdeqinout(100,alp,l,imu,enu,vm,bm,z,rofe,kc,nre,a,b,g0,gfn,
     .        psi4o(1,1,alp),psi4o(1,2,alp),psi4o(1,3,alp),psi4o(1,4,alp),
     .        psi4i(1,1,alp),psi4i(1,2,alp),psi4i(1,3,alp),psi4i(1,4,alp),valslo)
          enddo

C          call rdeqvs(1,l,imu,kc,nre,Enu,rofe,z,vm,bm,psi4o,sec,valslo(:,1,:))
C          call rdeqvs(1,l,imu,kc,nre,Enu,rofe,z,vm,bm,psi4i,sec,valslo(:,2,:))

C       (g2o - g1i - g2i) (Xo) = -g1o   M1
C       (f2o - f1i - f2i) (A)  = -f1o   M2
C       (G2o - G1i - G2i) (Xi) = -G1o   M3
C       (F2o - F1i - F2i)      = -F1o   M4

C         lam=1 : Satisfy M1,M3,M4, reserving M2 for energy search
          if (lam == 1) then
            sec(1,1) =  psi4o(kc,1,2)   !  g2o
            sec(2,1) =  psi4o(kc,3,2)   !  G2o
            sec(3,1) =  psi4o(kc,4,2)   !  F2o
            sec(1,2) = -psi4i(kc,1,1)   ! -g1i
            sec(2,2) = -psi4i(kc,3,1)   ! -G1i
            sec(3,2) = -psi4i(kc,4,1)   ! -F1i
            sec(1,3) = -psi4i(kc,1,2)   ! -g2i
            sec(2,3) = -psi4i(kc,3,2)   ! -G2i
            sec(3,3) = -psi4i(kc,4,2)   ! -F2i
            rhs(1) = -psi4o(kc,1,1)     ! -g1o
            rhs(2) = -psi4o(kc,3,1)     ! -G1o
            rhs(3) = -psi4o(kc,4,1)     ! -F1o

C         lam=2 : Satisfy M1,M2,M3, reserving M4 for energy search
          else
            sec(1,1) =  psi4o(kc,1,1)   !  g1o
            sec(2,1) =  psi4o(kc,2,1)   !  f1o
            sec(3,1) =  psi4o(kc,3,1)   !  G1o
            sec(1,2) = -psi4i(kc,1,2)   ! -g2i
            sec(2,2) = -psi4i(kc,2,2)   ! -f2i
            sec(3,2) = -psi4i(kc,3,2)   ! -G2i
            sec(1,3) = -psi4i(kc,1,1)   ! -g1i
            sec(2,3) = -psi4i(kc,2,1)   ! -f1i
            sec(3,3) = -psi4i(kc,3,1)   ! -G1i
            rhs(1)   = -psi4o(kc,1,2)   ! -g2o
            rhs(2)   = -psi4o(kc,2,2)   ! -f2o
            rhs(3)   = -psi4o(kc,3,2)   ! -G2o
          endif

          call dinv33(sec,0,seci,errmx)
          call dmpy31(0,seci,rhs,rhs)


C         psi4o(lam) -> psi4o(lam) + Xo psi4o(lamb)
          call daxpy(nre*4,rhs(1),psi4o(1,1,lamb),1,psi4o(1,1,lam),1)
C         psi4i(lam) -> A psi4i(lam) + Xi psi4i(lamb)
          call dscal(nre*4,rhs(2),psi4i(1,1,lam),1)
          call daxpy(nre*4,rhs(3),psi4i(1,1,lamb),1,psi4i(1,1,lam),1)
C         Components 1,3,4 should match to machine precision

C         Estimate correction to enu
          call rdeqvs(1,l,imu,kc,nre,enu,rofe,z,vm,bm,psi4o(1,1,lam),sec,valslo(:,1,:))
          call rdeqvs(1,l,imu,kc,nre,enu,rofe,z,vm,bm,psi4i(1,1,lam),sec,valslo(:,2,:))
C         Independent variable = enu
          Pv(1,1) = enu
C         Dependent variable = discontinuity in log derivative at kc
          Pv(2,1) = valslo(lam,1,2)/valslo(lam,1,1) - valslo(lam,2,2)/valslo(lam,2,1)
C         Estimate for de
          de = -valslo(lam,1,1)*valslo(lam,2,2) + valslo(lam,2,1)*valslo(lam,1,2)
          Pv(1:2,2) = Pv(1:2,1) ! Hang onto input P for printout

C          if (debug) then
C            call info8(1,0,0,
C     .      ' enu=%;9F  v(kc)=%;8F  s/v(o)=%;8F  s/v(i)=%;8F  de=%1,3;3g  delta(s/v)=%1,3;3g',
C     .      enu,valslo(lam,1,1),valslo(lam,1,2)/valslo(lam,1,1),
C     .        valslo(lam,2,2)/valslo(lam,2,1),de,Pv(2,1),7,8)
C
CC          call snit(iter,debug,gfn(:,alp,1,alp))
CC          if (iter < 0) cycle
C          endif

          if (ir /= 0 .or. iter == 1) then
C            if (debug) then
C              call info8(2,0,-1,
C     .        ' call rfalsi P %;18F   F %;10F  P-x1 %;10F P-x2 %;10F F1 %;10F F2 %;10F (P-x1)*(P-x2)<0=%l f1*f2<0=%l',
C     .        Pv,Pv(2,1),Pv(1,1)-wkg(1),Pv(1,1)-wkg(2),wkg(4),wkg(5),
C     .        (Pv(1,1)-wkg(1))*(Pv(1,1)-wkg(2)) < 0,wkg(4)*wkg(5) < 0)
C            endif
            call pshpr(ipr-70)
            call rfalsi(Pv,Pv(2,1),tolg,tolg,tolg,.03d0,30,wkg(1),ir)
            if (ir == 1) ir = 0
            if (ir == 0) norfalsi = .true.
            call poppr
C            if (debug) then
C              call info2(2,0,0,' New P bracket=%l ir=%i',(Pv(1,1)-wkg(1))*(Pv(1,1)-wkg(2)) < 0,ir)
C            endif
          endif

          userfalsi = (userfalsi .or. ir == -2 .or. ir == -3)
          if (norfalsi) userfalsi = .false.

C         New estimate for abscissa
          if (userfalsi) then
            enu = Pv(1,1)
          else
            enu = enu + de
          endif
          k1 = ir ; if (.not. userfalsi) k1 = -999; jpr = IPRT3; if (iter >= maxit-5) jpr = 10
          call info8(jpr,0,0,
     .      ' rdeqc  iter %i  %?#(n)#E#D# %;9F  deltaD %;9F  dE %;9F'//
     .      '  rfalsi=%?#n==-999#F #%-1j%i#  cnvg=%l  scl=%;3g',
     .      iter,1,Pv(1,2),Pv(2,2),de,k1,abs(de) < tole*enu,growth)
          if (abs(de) < abs(tole*enu)) exit

          if (iter == maxit) then
            call info8(10,1,0,
     .        ' rdeq00 exceeded maxit for l=%i, imu=%i, mu=%d, alp=%i  kc=%i  de=%d',
     .        l,imu,mu,alp,kc,de,7,8)
            call prrmsh('psi',rofe,psi4i(1,1,lam),nre,nre,4)
            call rx('rdeq failed to converge')
          endif

        enddo
        call info8(IPRT1,0,0,
     .    ' rdeqc cnvg l=%i, kap=%i,%25pmu=%d %33plam=%i  E=%;11F  nit=%i  d(slo/val)=%;9F  |dE|=%1;2g',
     .      l,kap,mu,lam,enu,iter,Pv(2,1),abs(de))

C   ... Assemble normalized psi(1:nre)
C       if (debug) then
C       call info2(10,0,0,' rdeq l imu %2,2i  psio(kc) %4;18,12D',[l,imu],psi4o(kc,:,lam))
C       call info2(10,0,0,' rdeq l imu %2,2i  psii(kc) %4;18,12D',[l,imu],psi4i(kc,:,lam))
C       endif
        do  k2 = 1, 4
          call dpcopy(psi4o(1,k2,lam),psi4i(1,k2,lam),1,kc,1d0)
        enddo
C       Normalize with quadrature from mesh containing all points 1:nr
        de = diracp(nre,psi4i(1,1,lam),psi4i(1,1,lam),rofi(1,2))
        call dscal(4*nre,1/dsqrt(de),psi4i(1,1,lam),1)

        icore = icore+1
        Ekmu(icore) = enu
C        if (debug) then
C          print *, 'imu,kap,mu,lam',imu,kap,mu,lam
C          call prrmsh('psi4i',rofe,psi4i(1,1,lam),nre,nre,4)
C        endif

C       stop 'not yet'
        call xyrhfr(Ekmu(icore),l,imu,z,nre,nre,psi4i(1,1,lam),rofe,rofe(1,2),vm,bm,rhom,amom,vrhom,tcorm,rhormxm)
        tcorD = tcorD + tcorm
        rhormxD = rhormxD + rhormxm
        if (opt1 >= 1) then
          do  i = 1, 2
            call daxpy(nre,1d0/2,rhom,1,rhoD(1,i),1)
            call daxpy(nre,dble(2*i-3)/2,amom,1,rhoD(1,i),1)
          enddo
        endif

        enddo                   ! Loop over lam
      endif                     ! 4 component case

      enddo                     ! Loop over imu
      call poppr

      if (icore /= 2*(2*l+1)) call rx('rdeqcore: core level count mismatch')

C     print *, 'integral of rho : ', ddot(nr,rhoD,1,rofi(1,2),1),ddot(nr,rhoD(1,2),1,rofi(1,2),1)

      deallocate(psi4o,psi4i,rofe,rhom,amom)

      end subroutine rdeqcore

      subroutine xyrhfr(ecore,l,imu,z,nr,nre,psi,rofi,wt,v,b,rho,amag,vrho,tcore,rhormx)
C- Make Dirac density and integrate potential*density for one core state
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
Ci   v     :Scalar part of spherical potential
Ci   b     :magnetic field.  Not used now
Co Outputs
Co   rho   :contribution core density from this state from 1..min(nr,nre)
Co   vrho  :integral
Co   rhormx:contribution to true density at nr from this state
Cr Remarks
Cr   xyrhfr makes the density at points 1..min(nr,nre), and the integral
Cr   of rho*density from 1..nre.  Thus:
Cr     if nre == nr, the density and integral are the same
Cr     if nre < nr  g**2 is negligible at nr (deep core state)
Cr     if nre > nr  (large sphere) contribution to rho*v is included
Cu Updates
Cu   10 Apr 12 Repackaged radial mesh quadrature (radwgt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nre,l,imu
      double precision z,ecore,vrho,tcore,rhormx
      double precision psi(nr,2,2),rofi(nr),wt(nr),v(nr),b(nr),rho(nr),amag(nr)
C ... Local parameters
      integer nrmx,ir,kap1,kap2
      double precision fpi,r,rhoir,rmax,momir
      double precision Rgg(2,2),Rff(2,2),mu,fac
C     double precision mqg1,mqg2,mqf1,mqf2
      procedure(real(8)) :: ddot

      call dpzero(rho,nr); call dpzero(amag,nr)

      fpi = 16d0*datan(1d0)
      vrho = 0
      nrmx = min0(nr,nre)
      mu  = imu - l - 1.5d0
      fac = 2d0*mu/(2*l+1d0)

C ... Make rho, and integrate vrho for points 1..nrmx
      do  ir = 2, nrmx
        r = rofi(ir)
        do  kap1 = 1, 2
          do  kap2 = 1, 2
            Rgg(kap1,kap2) = psi(ir,1,kap1)*psi(ir,1,kap2)
            Rff(kap1,kap2) = psi(ir,2,kap1)*psi(ir,2,kap2)
          enddo                 ! kap2
        enddo                   ! kap1
        rhoir = Rgg(1,1) + Rff(1,1) + Rgg(2,2) + Rff(2,2)
        momir = -fac*Rgg(1,1) + fac*Rgg(2,2) - 2d0*mu/(-2*l+1d0)*Rff(1,1) - 2d0*mu/(2*l+3d0)*Rff(2,2)

C        if ((imu == 1 .or. imu == 2*l+2) .and. .false.) then ! for a max imu (1 solution)
C          Rgg(2,2) = psi(ir,1,2)**2
C          Rff(2,2) = psi(ir,2,2)**2
C          rhoir = Rgg(2,2) + Rff(2,2)
C          momir = - 2d0*mu/(-2*l-1d0)*Rgg(2,2) + 2d0*mu/(2*l+3d0)*Rff(2,2)
C        else
C          do  kap1 = 1, 2
C            do  kap2 = 1, 2
C              Rgg(kap1,kap2) = psi(ir,1,kap1)*psi(ir,1,kap2)
C              Rff(kap1,kap2) = psi(ir,2,kap1)*psi(ir,2,kap2)
C            enddo               ! kap2
C          enddo                 ! kap1
C          rhoir =  Rgg(1,1) + Rff(1,1) + Rgg(2,2) + Rff(2,2)
C          mqg1 = -fac*Rgg(1,1) - 2d0*mu/(-2*l-1d0)*Rgg(2,2)
C          mqg2 = -dsqrt((2*l+1d0)**2-4*mu**2)/(2*l+1d0)*(Rgg(1,2)+Rgg(2,1))
C          mqf1 = -2d0*mu/(-2*l+1d0)*Rff(1,1) - 2d0*mu/(2*l+3d0)*Rff(2,2)
C          mqf2 = 0
C          momir = mqg1 + mqg2 - mqf1 - mqf2
C        endif
        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
        rho(ir) = rho(ir) + rhoir
        amag(ir) = amag(ir) + momir
      enddo
      rmax = rofi(nr)
      rhormx = 0
      if (nre >= nr) rhormx = rhoir/(fpi*rmax**2)

C ... Integrate rho*v from nrmx+1 .. nre
      do  ir = nrmx+1, nre
        call rx('xyrhfr not ready for nrmx>nre')
C        r = rofi(ir)
C        rhoir = g(ir,1)**2 + g(ir,2)**2 + g(ir,3)**2 + g(ir,4)**2
C        vrho = vrho + wt(ir) * rhoir * (v(ir)-2*z/r)
      enddo

      tcore = ecore - vrho

C     print *, 'xyrhfr : ', ddot(nre,rho,1,wt,1), ddot(nre,amag,1,wt,1)

      end
      subroutine rdeq0(E,z,vm,bm,rofi,nr,kc,a,b,l,imu,iopt,g0,gfn,psi,gsmt)
C- Solves radial Dirac equations for a given imu,l
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   vm    :spin-averaged spherical potential, extrapolated 2 radial mesh points
Ci   bm    :Magnetic field, extrapolated 2 radial mesh points
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   kc    :rofi(kc) is point where inwards and outwards solutions are matched
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l     :l quantum number
Ci   imu   :projection of j (positive integer, counts -l-1/2...l+1/2)
Ci   iopt  :1s digit = idmod
Ci         :0 linearization energy (input).  The full Dirac 4-vector is calculated.
Ci         :1 input E is an estimate only; its value is determined internally
Ci         :  as a function of the logarithmic derivative of the dominant component,
Ci         :  specified through gfn, only the log derivative of the dominant
Ci         :  large component is matched, providing an incomplete solution of the
Ci         :  Dirac equation.
Ci         :2 similar to iopt=1, except that the cross-coupling between (g1,f1) and
Ci         :  (g2,f2) through the magnetic field is omitted.  This ensures that
Ci         :  (g2,f2) vanish for alpha=1 and (g1,f1) vanish for alpha=2; the result
Ci         :  satisfies an approximate differential equation
Ci   g0    :Initial guess for value of dominant large component g_alpha at first point.
Ci         :If 0, value is chosen internally.
Cio Inputs/Outputs
Cio  E     :eigenenergy
Cio        :If idmod=0, enu is input, and psi is calculated for given E
Cio        :All four components are determined for each pair of solns alpha=1,2
Cio        :If idmod=1, input enu is an estimate only;
Cio        :slope/value at rofi(nr) is fixed.  See description of idmod in opt, above.
Cio        :Output E depends alpha; E(1,2) are returned for each alpha
Cio  gfn   :Boundary conditions for Dirac equation, 8 components to include alpha=1,2
Cio                         value            energy derivative
Cio          g1(nr)     gfn(1,1,1,alpha)     gfn(1,1,2,alpha)
Cio          f1(nr)     gfn(2,1,1,alpha)     gfn(2,1,2,alpha)
Cio          g2(nr)     gfn(1,2,1,alpha)     gfn(1,2,2,alpha)
Cio          f2(nr)     gfn(2,2,1,alpha)     gfn(2,2,2,alpha)
Cio        :gfn(1:2,1:2,1,alpha) are the boundary conditions; only 4 used for a given alpha
Cio        :gfn(1:2,1:2,2,alpha) are approximate energy derivatives, and are not touched
Cio        :                     They are used to facilitate convergence
Cio        :On input, if zero, an estimate for gfn is taken from the scalar Dirac equation
Cio        :On input, if nonzero, gfn is used as a starting guess for the boundary conditions
Cio        :On output,  boundary conditions for the given v and enu, for the given l,imu
Co Outputs
Co   psi   :psi(meshpoint,1,i,alpha) large component g
Co         :psi(meshpoint,2,i,alpha) small component f
Co         :i = 1 or 2 corresponding to kappa or kappa' = -kappa-1.
Co         :alpha = 1 or 2, called lambda in Schick's paper (PhysRevB, 54 1610 (1996))
Co   gsmt  :gsmt(i,alpha) = slope dg/dr of g_i,alpha at augmentation radius, i=1,2 alpha=1,2
Cl Local variables
Cl   valslo: values and slopes of large components gi at matching point kc
Cl         : valslo(i,j,k) : i=1,2 labels kappa; j=1,2 for (in,out); k=1,2 for (val,slo)
Cl         : valslo(i,1,1) valslo(i,2,1)  valslo(i,1,2) valslo(i,2,2)
Cl         :    g_i(out)     g_i(in)       g'_i(out)     g'_i(in)
Cr Remarks
Cr   gn(alpha,1:2,) and fn(alpha,1:2) must be solved as four coupled equations
Cr   except when imu=1 or imu=2l*2.  Then gn(alpha,2-alpha) and fn(alpha,2-alpha) vanish.
Cr
Cr   In this version enu depends only on kappa
Cu Updates
Cu   01 Aug 16 (MvS) Redesign of solver
Cu   30 Jul 16 (MvS) First attempt at redesign of solver
Cu   24 Oct 14 (Vishina) Rewrite of rseq.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,imu,l,kc,iopt
      double precision z,a,b,E(2),g0(2)
      double precision psi(nr,2,2,2),vm(nr+2),bm(nr+2),rofi(nr,2)
      double precision gsmt(2,2),gfn(2,2,2,2)
C ... Local parameters
C     logical :: debug = .false.
      logical lcnvg,lcnvg2,userfalsi,getbc,norfalsi
      integer stdo,ir,iter,alp,ii,jj,ipr,ipath,idmod,modrio
      double precision enu,c,mu,normf,ampg,ampf,csq,de,gam,q,growth
      double precision pvsav(4),wkg(28),gfx(2,2,2,2),gfnl(2,2,2,2)
      double precision kap,kapp,rmax,absg,absdP,valslo(2,2,2),vsnr(2,2)
      real(8) :: Pv(4,10),Jac(4,4)
      real(8) :: psi4o(nr,4),psi4i(nr,4),psi4ox(nr,4),psi4ix(nr,4)
      real(8), parameter :: dele=.006d0  ! For finite difference making psidot
      real(8), parameter :: tole=1d-12   ! Tolerance in energy
      real(8), parameter :: tolg=1d-13   ! Tolerance in w.f. discontinuity
      real(8), parameter :: tolg2=1d-11  ! Softer Tolerance in w.f. discontinuity
      real(8), parameter :: grfac = 1.1d0 ! Growth factor in extrapolating to root
      integer, parameter :: IPRT1=60, IPRT2=70, IPRT3=100, maxit=100
C     integer, parameter :: IPRT1=60/10, IPRT2=100/10, IPRT3=100, maxit=200*1
      procedure(integer) :: idamax,nglob,iprint
      procedure(real(8)) :: dot3,dlength,diracp

C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      data wkg /28*0d0/

C ... Setup
      if (z == 0) return
      enu = E(1)
      gfnl = gfn                ! Preserve local copy of initial b.c. for debugging
      idmod = 0; if (iopt /= 0) idmod = 1 ! Whether to freeze enu or slo/val
      modrio = 100; if (iopt == 2) modrio = 101 ! call rdeqinout in this mode
C     savenu = enu
      call getpr(ipr)
      call pshpr(ipr)           ! Preserve so we can locally modify it
      stdo = nglob('stdo')
      mu  = imu - l - 1.5d0
      rmax = rofi(nr,1)
      csq = c*c
      kap = l; kapp = -l-1; if (imu == 1 .or. imu == 2*l+2) kap = -l-1


C --- Big loop over independent solutions alp ---
      getbc = dlength(size(gfn),gfn,1) == 0
      do  alp = 1, 2

C     if (l == 2 .and. imu == 1) call snot(1,debug)
C     if (l == 2 .and. imu == 2 .and. alp == 2) call snot(1,debug)
C     if (l == 1 .and. imu == 2) call snot(1,debug)

C ... Integrate outward to obtain initial boundary conditions for inwards integration
      ii = 2*alp-1
C     if (idmod == 0) then
      if (getbc) then
        do  jj = 0, 1  ! 0, 1 for energy derivative
          call rdeqinout(10+modrio,alp,l,imu,enu+jj*dele,vm,bm,z,rofi,0,nr,a,b,g0,gfn,
     .      psi4o(1,1),psi4o(1,2),psi4o(1,3),psi4o(1,4),
     .      psi4i(1,1),psi4i(1,2),psi4i(1,3),psi4i(1,4),valslo)
C         if (jj==0) call prrmsh('psi',rofi,psi4o,nr,nr,4)
C         Approximately normalize
          call prodrw(nr,2,psi4o(1,ii),psi4o(1,ii),rofi(1,2),q); q=dsqrt(q)
          call dscal(nr*2,1/q,psi4o(1,ii),1)
          g0(alp) = psi4o(2,ii) ! g(2), approximately normalized
          gfn(1:2,alp,jj+1,alp) = psi4o(nr,ii:ii+1)
        enddo
        gfn(:,alp,2,alp) = (gfn(:,alp,2,alp) - gfn(:,alp,1,alp))/dele
      endif

      call rdeqvs(0,l,imu,nr,nr,enu,rofi,z,vm,bm,psi4i,gfn(1,1,1,alp),vsnr)
      gam = dsqrt(kap**2 - (2*Z/c)**2)
      if (alp == 2) gam = dsqrt(kapp**2 - (2*Z/c)**2)
      call info8(IPRT2,1,-1,
     .  ' rdeq0 start l=%i, imu=%i, mu=%d, alp=%i, idmod=%i  E=%;9F  D=%;9F  gamma=%1,3;3g',
     .    l,imu,mu,alp,idmod,enu,rmax*vsnr(alp,2)/vsnr(alp,1),gam)
      call info2(IPRT2,0,0,'  nr=%i',nr,2)
C ... Construct wave function satisfying bc only for (g1,f1), alp=1 or for (g2,f2), alp=2
C ... 1 normalize initial guess at psi;
C     2 Make dominant (f,g) ratio consistent with energy, neglecting other boundary cond.
C       If idmod=1, f(nr)/g(nr) is fixed.  Otherwise enu is fixed.
      iter = 0; ir = 0; userfalsi = .false.; growth = 1; ipath = 0; ampf = 1
      norfalsi = .false.  ! When true, rfalsi acceleration is not used
      do  while (iter <= maxit)
        iter = iter+1
        ii = 2*alp-1
        if (iter >= 30) norfalsi = .true.
        if (iter >= maxit-5) then
          call setpr(max(iprint(),IPRT3))
          ipr = iprint()
        endif

C        if (debug) then
C          call rdeqvs(0,l,imu,nr,nr,enu,rofi,z,vm,bm,psi4i,gfn(1,1,1,alp),vsnr)
C          call info5(1,0,0,' call rdeq g,f,ratio %2:2;9F  %;9F  v,s,D %2:2;9F  %;9F',
C     .      gfn(1,alp,1,alp),gfn(2,alp,1,alp)/gfn(1,alp,1,alp),
C     .      vsnr(alp,:),rmax*vsnr(alp,2)/vsnr(alp,1),5)
C        endif

        call rdeqinout(modrio,alp,l,imu,enu,vm,bm,z,rofi,kc,nr,a,b,g0,gfn,
     .    psi4o(1,1),psi4o(1,2),psi4o(1,3),psi4o(1,4),
     .    psi4i(1,1),psi4i(1,2),psi4i(1,3),psi4i(1,4),valslo)

C       Join dominant gin,gout at kc to make single function; normalize approximately
        normf = valslo(alp,2,1)/valslo(alp,1,1)
C       Scale by -1 if psi4i*psi4o < 0
        if (normf < 0) then
          call info5(IPRT1,0,0,
     .    ' rdeq0 (warning) l=%i, imu=%i, mu=%d, alp=%i normf changed sign',l,imu,mu,alp,5)
C          print *, kc, nr, rofi(kc,1), g0(alp)
C          call prrmsh('psio',rofi,psi4o,nr,nr,4)
C          call prrmsh('psii',rofi,psi4i,nr,nr,4)
          call dscal(2,-1d0,gfn(1,alp,1,alp),1)
          call dscal(2,-1d0,valslo(alp,2,1:2),1)
          call dscal(4*nr,-1d0,psi4i,1)
          normf = -normf
          norfalsi = .true.
        endif
C       Scale psi4o to match corresponding large component of psi4i
C       call dpcopy(psi4o(1,ii),psi4i(1,ii),1,kc,normf)
C       call dpcopy(psi4o(1,ii+1),psi4i(1,ii+1),1,kc,normf)
        do  jj = 1, 4
          call dpcopy(psi4o(1,jj),psi4i(1,jj),1,kc,normf)
        enddo
        valslo(alp,1,:) = valslo(alp,1,:) * normf ! scale v+s(out) to
C       Normalize
        call prodrw(nr,2,psi4i(1,ii),psi4i(1,ii),rofi(1,2),q); q=dsqrt(q)
        valslo(alp,:,:) = valslo(alp,:,:)/q
        gfn(:,alp,1,alp) = gfn(:,alp,1,alp)/q
        call dscal(nr*2,1/q,psi4i(1,ii),1)
        g0(alp) = psi4i(2,ii)     ! g(2), approximately normalized

C   ... Extract abscissa, ordinate for given boundary conditions
C       Independent variable = log derivative of g at kc
        call rdeqvs(1,l,imu,nr,nr,enu,rofi,z,vm,bm,psi4i,gfn(1,1,1,alp),vsnr)
        Pv(1,1) = rmax*vsnr(alp,2)/vsnr(alp,1)
        if (idmod == 1) Pv(1,1) = enu
C       Dependent variable = discontinuity in log derivative at kc
        Pv(2,1) = valslo(alp,1,2)/valslo(alp,1,1) - valslo(alp,2,2)/valslo(alp,2,1)
        if (iter == 1) Pv(1:2,2) = Pv(1:2,1)  ! No prior pair yet ... probably not needed
        if (iter == 1) Pv(1:2,3) = Pv(1:2,1)  ! Keep starting pair ... probably not needed
        Pv(1:2,4) = Pv(1:2,1)                 ! Hang onto input P for printout

C        if (debug) then
C        call info2(2,0,0,' P now   %;10F   F now   %;10F',Pv,Pv(2,1))
C        endif

C        if (debug) then
C          call info8(1,0,0,
C     .      ' Pv=%2:2;9F  g0=%;8F  normf=%;8F  q=%;8F  v(kc)=%;8F  s/v(o)=%;8F  s/v(i)=%;8F gfn=%2:2;9F',
C     .      Pv,g0(alp),normf,q,valslo(alp,1,1),valslo(alp,1,2)/valslo(alp,1,1),
C     .      valslo(alp,2,2)/valslo(alp,2,1),gfn(:,alp,1,alp))
CC          call snit(iter,debug,gfn(:,alp,1,alp))
CC          if (iter < 0) cycle
C        endif

C       Estimate energy shift de = g(in)*(g'(in) - g'(out)), if g is normalized
        de = -valslo(alp,2,1)*(valslo(alp,2,2) - valslo(alp,1,2))

C       Get rfalsi's estimate for new abscissa D
        if (ir /= 0 .or. iter == 1) then
C        if (debug) then
C          call info8(2,0,-1,
C     .    ' call rfalsi P %;18F   F %;10F  P-x1 %;10F P-x2 %;10F F1 %;10F F2 %;10F (P-x1)*(P-x2)<0=%l f1*f2<0=%l',
C     .    Pv,Pv(2,1),Pv(1,1)-wkg(1),Pv(1,1)-wkg(2),wkg(4),wkg(5),
C     .    (Pv(1,1)-wkg(1))*(Pv(1,1)-wkg(2)) < 0,wkg(4)*wkg(5) < 0)
C        endif
        call pshpr(ipr-70)
        call rfalsi(Pv,Pv(2,1),tolg,tolg,tolg,.03d0,30,wkg(1),ir)
        if (ir == 1) ir = 0
        if (ir == 0) norfalsi = .true.
        call poppr
C        if (debug) then
C          call info2(2,0,0,' New P bracket=%l ir=%i',(Pv(1,1)-wkg(1))*(Pv(1,1)-wkg(2)) < 0,ir)
C        endif
        endif

        userfalsi = (userfalsi .or. ir == -2 .or. ir == -3)
!       if (abs(de) < tole) norfalsi = .true.
        if (norfalsi) userfalsi = .false.

C   ... New estimate for energy (independent variable)
        if (idmod == 1) then
          if (userfalsi) then
            enu = Pv(1,1)
          else
            enu = enu + de
          endif
C   ... Estimate shift in g, f to keep de=0
        else
          if (userfalsi) then
            vsnr(alp,1) = gfn(1,alp,1,alp)
            vsnr(alp,2) = vsnr(alp,1)*Pv(1,1)/rmax  ! Translate new abscissa into b.c. at nr
            call rdeqvs(10,l,imu,nr,nr,enu,rofi,z,vm,bm,psi4i,gfn(1,1,1,alp),vsnr)
          else
            if (iter > 70) ipath = 1
C           Normal case : vary both g and f at augmentation radius
            if (ipath == 0) then
              gfn(:,alp,1,alp) = gfn(:,alp,1,alp) - growth * de * gfn(:,alp,2,alp)
C           Safe case: vary only f at augmentation radius
            else
              gfn(2,alp,1,alp) = gfn(2,alp,1,alp) - growth * de * gfn(2,alp,2,alp)
            endif
            if (Pv(2,1)*Pv(2,2) > 0 .and. iter > 1) growth = min(growth*grfac,2d0)
            if (Pv(2,1)*Pv(2,2) <= 0) growth = 1
C           print *, 'iter, growth, q', iter,growth,q
          endif
        endif
        jj = ir ; if (.not. userfalsi) jj = -999
        call info8(IPRT3,0,0,
     .    ' rdeq00  iter %i  %?#(n)#E#D# %;9F  deltaD %;9F  dE %;9F'//
     .    '  rfalsi=%?#n==-999#F #%-1j%i#  cnvg=%l  scl=%;3g',
     .    iter,idmod,Pv(1,4),Pv(2,4),de,jj,abs(de) < tole,growth)
        if (abs(de) < tole) exit
        if (iter >= maxit-5 .and. abs(de/e(alp)*10) < tole) exit

        Pv(:,2) = Pv(:,1)       ! Current pair becomes prior pair

C        if (debug) then
C          call snit(iter,debug,gfn(:,alp,1,alp))
C          if (iter < 0) cycle
C        endif

      enddo

      if (iter > maxit) then
        call info8(10,1,0,
     .    ' rdeq00 exceeded maxit for l=%i, imu=%i, mu=%d, alp=%i  kc=%i  de=%d',
     .    l,imu,mu,alp,kc,de,7,8)
        call prrmsh('psi',rofi,psi4i,nr,nr,4)
        call rx('rdeq failed to converge')
      endif
      call info8(IPRT1,0,0,
     .  ' rdeq00 cnvg l=%i, imu=%i, mu=%d, alp=%i  nit %,2i   E=%;9F  D=%;9F  |dE|=%1;2g',
     .    l,imu,mu,alp,iter,enu,rmax*vsnr(alp,2)/vsnr(alp,1),abs(de))

C ... Early exit when no coupling between kap1 and kap2
      if (imu == 1 .or. imu == 2*l+2 .or. iopt > 0) then
        call dcopy(4*nr,psi4i,1,psi(1,1,1,alp),1)
        if (idmod == 1) E(alp) = enu
        if (imu == 1 .or. imu == 2*l+2) then  ! alp=1 and alp=2 are equivalent
          call dpzero(psi(1,1,1,2),nr*2)
          call dcopy(nr*2,psi(1,1,1,1),1,psi(1,1,2,2),1)
          exit
        endif
        cycle   ! When alp=1 and alp=2 solutions may be different
      endif

C     if (debug) call pshpr(70)

C ... Integrate outward to obtain initial boundary conditions for inwards integration
      call rdeqinout(10+modrio,alp,l,imu,enu,vm,bm,z,rofi,0,nr,a,b,g0,gfn,
     .  psi4o(1,1),psi4o(1,2),psi4o(1,3),psi4o(1,4),
     .  psi4i(1,1),psi4i(1,2),psi4i(1,3),psi4i(1,4),valslo)

      if (alp == 1) then
C       gfn(1,1,1,1) = psi4o(nr,1) ! Use bc derived from 2 component solution
C       gfn(2,1,1,1) = psi4o(nr,2) ! Use bc derived from 2 component solution
        gfn(1,2,1,1) = psi4o(nr,3)
        gfn(2,2,1,1) = psi4o(nr,4)
      else
        gfn(1,1,1,2) = psi4o(nr,1)
        gfn(2,1,1,2) = psi4o(nr,2)
C       gfn(1,2,1,2) = psi4o(nr,3) ! Use bc derived from 2 component solution
C       gfn(2,2,1,2) = psi4o(nr,4) ! Use bc derived from 2 component solution
      endif

C --- Loop to find for F[P] which gives F[P]=0.  F is stored in Pv(:,2) ---
      ir = 0; iter = 0; pvsav = 0; ampg = 1; ampf = 1
      do  while (iter <= maxit)
      iter = iter+1

C ... Dirac equation with trial boundary conditions P
      ii = kc ; if (iter == 1) ii = 0
      call rdeqinout(modrio,alp,l,imu,enu,vm,bm,z,rofi,ii,nr,a,b,g0,gfn,
     .  psi4o(1,1),psi4o(1,2),psi4o(1,3),psi4o(1,4),
     .  psi4i(1,1),psi4i(1,2),psi4i(1,3),psi4i(1,4),valslo)

C ... Make F[P]
C     Initial value for P : gfn scaled by ampg (g components) and ampf (f components)
      if (iter == 1) then
        ampg = 1/max(abs(psi4o(nr,1)),abs(psi4o(nr,3)))
        ampf = 1/max(abs(psi4o(nr,2)),abs(psi4o(nr,4)))
        call gftopv(.true.,imu,l,alp,ampf,ampg,gfn,Pv) ! Map gfn to Pv
      endif
      Pv(1,2) = (psi4o(kc,1) - psi4i(kc,1))*ampg
      Pv(2,2) = (psi4o(kc,2) - psi4i(kc,2))*ampf
      Pv(3,2) = (psi4o(kc,3) - psi4i(kc,3))*ampg
      Pv(4,2) = (psi4o(kc,4) - psi4i(kc,4))*ampf

C ... Re-entry point for new line minimization (ir=-2, in which case P is unchanged from last call)
   12 continue

      absg = dlength(4,Pv(1,2),1) ! |F|
      call info2(IPRT3,1,0,' P now %4:2;12F  iter %i',Pv,iter)
      call info2(IPRT3,0,0,' F now %4:2;12F  |F|=%1;3g',Pv(1,2),absg)
      pvsav(1:4) = Pv(1:4,1) ! Pv may be updated with
      lcnvg = absg < tolg  ! Convergence criterion
      if (lcnvg) exit

C ... Case new line minimization : make Jacobian Jij = dpsi_j(kc)/dgfn(i)
      if (ir == 0 .or. ir == -2) then
        do  jj = 1, 4
          Pv(jj,1) = pvsav(jj) + 1/2d0**20
          call gftopv(.false.,imu,l,alp,ampf,ampg,gfx,Pv) ! Map Pv to temporary gfn
          call rdeqinout(modrio,alp,l,imu,enu,vm,bm,z,rofi,kc,nr,a,b,g0,gfx,
     .      psi4ox(1,1),psi4ox(1,2),psi4ox(1,3),psi4ox(1,4),
     .      psi4ix(1,1),psi4ix(1,2),psi4ix(1,3),psi4ix(1,4),valslo)
          Jac(1,jj) = (psi4ox(kc,1) - psi4ix(kc,1))*ampg
          Jac(2,jj) = (psi4ox(kc,2) - psi4ix(kc,2))*ampf
          Jac(3,jj) = (psi4ox(kc,3) - psi4ix(kc,3))*ampg
          Jac(4,jj) = (psi4ox(kc,4) - psi4ix(kc,4))*ampf
          do  ii = 1, 4
            Jac(ii,jj) = (Jac(ii,jj) - Pv(ii,2))*2d0**20
          enddo
          Pv(jj,1) = pvsav(jj)
        enddo
      endif

C       call gradzr(4,Pv,Jac,1d-5,1d0/2,tolg,tolg,.3d0,wkg,320,ir)
      call pshpr(ipr-70)
      call gradzr(4,Pv,Jac,1d-5,1d0,tolg,tolg,.3d0,wkg,320,ir)
      call poppr
      if (ir > 0) call rx1('rdeq: gradzr returned ir=%i ... aborting',ir)
      absdP = dlength(4,pvsav(1:4)-Pv(1:4,1),1)
      lcnvg2 = absdP < tolg .and. absg < tolg2 ! Allow softer tolerance in g if change in P is also small
      if (lcnvg2 .and. ir /= -2) ir = 0

      call info8(IPRT3,0,0,' rdeq0 iter %i ln %d: imu=%i l=%i %1j'//
     .  '%-1j%?#n==-0#cnvg##'//
     .  '%-1j%?#n==-1#continue line##'//
     .  '%-1j%?#n==-2#new line##'//
     .  '%-1j%?#n>1#trouble!!##  dxmx=%;7F  |g|=%1;3g  |dP|=%1;3g',
     .    iter,wkg(19),imu,l,ir,wkg(1)*Pv(idamax(4,Pv(1,4),1),4),absg,absdP)

      call gftopv(.false.,imu,l,alp,ampf,ampg,gfn,Pv)  ! Map Pv to gfn

      if (ir == -2 .and. dlength(4,pvsav(1:4)-Pv(1:4,1),1) == 0) goto 12 ! Pv unchanged
      if (iter > maxit) call fexit3(-1,111,' Exit -1 rdeq0 l=%i, imu=%i:  '//
     .      'failed to converge to tolerance.  |g|=%1;3g',l,imu,absg)
      if (ir == 0) exit   ! Iterate until converged
      enddo

C ... Solution has converged
      if (dlength(4,pvsav(1:4)-Pv(1:4,1),1) /= 0) then
        call rdeqinout(modrio,alp,l,imu,enu,vm,bm,z,rofi,kc,nr,a,b,g0,gfn,
     .    psi4o(1,1),psi4o(1,2),psi4o(1,3),psi4o(1,4),
     .    psi4i(1,1),psi4i(1,2),psi4i(1,3),psi4i(1,4),valslo)
      endif

      ii = 2*alp-1
      call rdeqvs(1,l,imu,nr,nr,enu,rofi,z,vm,bm,psi4i,gfn,vsnr)
      call info8(IPRT1,0,0,
     .  ' rdeq0  cnvg l=%i, imu=%i, mu=%d, alp=%i  nit %,2i  E=%;9F  D=%;9F  |dF|=%1;2g',
     .  l,imu,mu,alp,iter,enu,rmax*vsnr(alp,2)/vsnr(alp,1),absg)

C ... Accumulate 4-vector from outward and inwards integration
      psi(1:kc,1,1,alp) = psi4o(1:kc,1) !g1, outward integration
      psi(1:kc,2,1,alp) = psi4o(1:kc,2) !f1, outward integration
      psi(1:kc,1,2,alp) = psi4o(1:kc,3) !g2, outward integration
      psi(1:kc,2,2,alp) = psi4o(1:kc,4) !f2, outward integration
      psi(kc+1:nr,1,1,alp) = psi4i(kc+1:nr,1) !g1, inward integration
      psi(kc+1:nr,2,1,alp) = psi4i(kc+1:nr,2) !f1, inward integration
      psi(kc+1:nr,1,2,alp) = psi4i(kc+1:nr,3) !g2, inward integration
      psi(kc+1:nr,2,2,alp) = psi4i(kc+1:nr,4) !f2, inward integration

      enddo ! Loop over alp

C ... Render psi(alpha=2) and psi(alpha=1) orthonormal
      if (imu /= 1 .and. imu /= 2*l+2) then
C       Orthogonalize psi(alpha=2) to psi(alpha=1); normalize both
C       call orthdirac(11,nr,rofi(1,2),psi,psi(1,1,1,2))
C       Symmetrically orthogonalize psi(alpha=1) and psi(alpha=2); normalize both
        call orthdirac(12,nr,rofi(1,2),psi,psi(1,1,1,2))
      endif

C     grep 'rdeq norm' out | awk '{print $10, $11}' | mc . -e1 x1-x2 -px
C      call info5(10,0,0,' rdeq norm l imu %2,2i  <ai aj> %3;14,9D  s1  s2  %;12F  %;12F',
C     .  [l,imu],
C     .  [diracp(nr,psi,psi,rofi(1,2)),
C     .   diracp(nr,psi,psi(1,1,1,2),rofi(1,2)),
C     .   diracp(nr,psi(1,1,1,2),psi(1,1,1,2),rofi(1,2))],
C     .   sum(psi(:,:,:,1)),sum(psi(:,:,:,2)),0)

C ... Ensure dominant parts of alpha=1 and alpha=2 solutions have the same phase at RMT.
      if ((psi(nr,1,1,1)*psi(nr,1,2,2)) < 0) then
        psi(:,:,:,2) = -psi(:,:,:,2)
      endif

C fplot -map 1:nr -colsy 2 out.ni -map 1:nr -ord 'x3*274' -lt 1,col=1,0,0 out.ni   -lt 2,col=0,0,0 -ord 'x4*274' out.ni   -lt 2,col=1,0,0 -ord 'x5*274^2' out.ni
C      if (l == 1 .and. imu == 2) then
C        print *, 'kc=',kc
C        call prrmsh('psi',rofi,psi,nr,nr,4)
C      endif

C ... For the potential parameters (maybe not needed ?)
C      kap = l; kapp = -l-1
C      if (imu == 1 .or. imu == 2*l+2) kap = -l-1
C      up = mu/(l-0.5d0) ; u = mu/(l+0.5d0) ; upp = mu/(l+1.5d0)
C      if (imu == 1 .or. imu == 2*l+2) then
C        nug1 = 1d0 + 2*z/rmax/csq + (enu - vm(nr) - upp*bm(nr))/csq
C      else
C        nug1 = 1d0 + 2*z/rmax/csq + (enu - vm(nr) + up *bm(nr))/csq
C      endif
C      nug2 = 1d0 + 2*z/rmax/csq + (enu - vm(nr) - upp*bm(nr))/csq
C      gsmt(1,:) = (c*nug1*psi(nr,2,1,:) - kap /rmax*psi(nr,1,1,:))
C      gsmt(2,:) = (c*nug2*psi(nr,2,2,:) - kapp/rmax*psi(nr,1,2,:))
      do  alp = 1, 2
        call rdeqvs(1,l,imu,nr,nr,enu,rofi,z,vm,bm,psi(1,1,1,alp),gfn,vsnr)
        gsmt(1:2,alp) = vsnr(1:2,2)
      enddo

C     if (debug) stop 'rdeq0'
C     if (debug) call poppr

      call info2(IPRT3,0,0,' slope of g at augmentation radius: %4:2;12F',gsmt,2)
      call poppr

      end subroutine rdeq0

      subroutine rdeqinout(mode,alpha,l,imu,enu,vm,bm,z,rofi,kc,nr,a,b,g0,gfn,
     .                     g1out,f1out,g2out,f2out,g1in,f1in,g2in,f2in,valslo)
C- Integrates radial Dirac equation both inwards and outwards
C ----------------------------------------------------------------------
Ci   mode  :1s digit
Ci         :1 Eliminate coupling linking (g1,f1) to (g2,f2) by setting sqrt(u)=0
Ci         :10s digit
Ci         :1 Outward integration only
Ci         :2 Inward integration only
Ci         :100s digit
Ci         :0  Do outwards integration numerically on psi
Ci         :1  Do outwards integration numerically on (r^-gamma psi)
Ci   alpha :1 or 2 - labels two possible independent solutions of Dirac equation.
Ci   l     :l quantum number
Ci   imu   :projection of j (positive integer, counts -l-1/2...l+1/2) (mu  = imu - l - 1.5d0)
Ci   enu   :Linearization energy
Ci   vm    :potential
Ci   bm    :magnetic field
Ci   z     :nuclear charge
Ci   rofi  :radial mesh points
Ci   kc    :final mesh point for inward and outward solutions, and where valslo is evaluated
Ci         :kc=0 => solutions generated for entire mesh, valslo evaluated at nr
Ci         :kc<0 => solutions generated for entire mesh, valslo evaluated at -kc
Ci   nr    :number of radial mesh points
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   g0    :value of g_alpha at first point.  If 0, value is chosen internally.
Ci         :f_alpha, g_(2-alpha) and f_(2-alpha) are determined from small-r
Ci         :behavior of Dirac equation
Ci   gfn   :Boundary conditions for Dirac equation at rofi(nr).
Ci         :Only four of them are used, depending on alpha.
Ci                          value
Ci           g1(nr)    gfn(1,1,1,alpha)
Ci           f1(nr)    gfn(2,1,1,alpha)
Ci           g2(nr)    gfn(1,2,1,alpha)
Ci           f2(nr)    gfn(2,2,1,alpha)
Co Outputs
Co   ... integrated from the origin outwards
Co   g1out : Large lambda=1 component of wave function times r.  Dominant for alpha=1
Co   g2out : Large lambda=2 component of wave function times r.  Dominant for alpha=2
Co   f1out : Small lambda=1 component of wave function times r
Co   f2out : Small lambda=2 component of wave function times r
Co   ... integrated from the MT boundary inwards
Co   g1in  : analog of g1out
Co   f1in  : analog of f1out
Co   g2in  : analog of g2out
Co   f2in  : analog of f2out
Co   valslo: values and slopes of large components gi at matching point kc
Co         : valslo(i,j,k) : i=1,2 labels kappa; j=1,2 for (in,out); k=1,2 for (val,slo)
Co         : valslo(i,1,1) valslo(i,2,1)  valslo(i,1,2) valslo(i,2,2)
Co         :    g_i(out)     g_i(in)       g'_i(out)     g'_i(in)
Cl Local variables
Cr Remarks
Cr    This routine integrates one of two independent solutions to the Dirac equation.
Co    g2, f2 are coupled to g1,f2 by magnetic field.  One pair is dominant;
Co    the other pair is small and vanishes if B=0, or if imu= 1 or 2*l+2.
Cr    Independent solutions alpha=1,2 differ through their boundary conditions.
Cr    For small r,
Cr      alpha=1 => g1 and f1 are dominant, proportional to r^gam1 for small r
Cr                 g2 and f2 go to zero at a rate r^(gam1+1), and vanish if B=0
Cr      alpha=2 => g2 and f2 are dominant, proportional to r^gam2 for small r
Cr                 g1 and f1 go to zero at a rate r^(gam2+1), and vanish if B=0
Cr    Here
Cr      gam1 = sqrt(kap1**2 - (2*Z/c)**2)         gam2 = sqrt(kap2**2 - (2*Z/c)**2)
Cr      kap1 = l, or -l-1 if imu=1 or imu=2l+2    kap2 = -l-1
Cr      kap  = l : spin and orbital antialigned;  kap  = -l-1 spin and orbital aligned
Cr    See rdinit below for more detailed description of small r behavior
Cr
Cr    Dirac equations are:
Cr      g1'(r) = -kap1/r g1(r) + nug1(r) (c f1(r))                     (D1)
Cr      f1'(r) = -kap1/r f1(r) + nuf1(r) g1(r) - B(r)/c sqrt(u) g2(r)  (D2)
Cr      g2'(r) = -kap2/r g2(r) + nug2(r) (c f2(r))                     (D3)
Cr      f2'(r) = -kap2/r f2(r) + nuf2(r) g2(r) - B(r)/c sqrt(u) g1(r)  (D4)
Cr    where
Cr      nug1(r) =  (1 + (e1 + 2Z/r - V(r) + ug1*B(r))/c^2
Cr      nuf1(r) = -(0 + (e1 + 2Z/r - V(r) + uf1*B(r))
Cr      nug2(r) =  (1 + (e2 + 2Z/r - V(r) + ug2*B(r))/c^2
Cr      nuf2(r) = -(0 + (e2 + 2Z/r - V(r) + uf2*B(r))
Cr    and
Cr      kap2 = kapp = -l-1
Cr      ug2  = -upp = -mu/(l+1.5d0)
Cr      uf2  = -u   = -mu/(l+0.5d0)
Cr       e2  = enu2
Cr    and if imu=1 or imu=2l+2
Cr      (kap1, ug1, uf1, e1) = (kap2, ug2, uf2, e2)
Cr    and if not
Cr      kap1 = l
Cr      ug1  = +up =   mu/(l-0.5d0)
Cr      uf1  = +u  =   mu/(l+0.5d0)
Cr       e1  = enu1
Cr    This routine assumes enu1=enu2=enu
Cr    Compare to Turek's book Eq. 6.13, 6.14 (case B(r)=0)
Cr      (d/dr + (1+k)/r) gtrue_k - (1+(E-V(r))/c^2) cftrue_k = 0
Cr      (d/dr + (1-k)/r) ftrue_k +   (E-V(r)) gtrue_k        = 0
Cr    Note (g,f) = r*(gtrue,ftrue)
Cr  Use : r * d/dr (gtrue_k) = d/dr g_k - g_k/r
Cr  Substituting (g,f) for (gtrue,ftrue):
Cr      (d/dr + k/r) g_k - (1+(E-V(r))/c^2) c*f_k = 0  which is (D1) when B(r)=0
Cr      (d/dr - k/r) f_k +    (E-V(r)) g_k        = 0  which is (D2) when B(r)=0
Cu Updates
Cu   20 Aug 16 Bug fix, inward integration
Cu   01 Aug 16 Redesign with better boundary conditions; integrate r^-gam * (g,f)
Cu   29 Sep 14 (Alena Vishina) rewrite of rseq.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,alpha,l,imu,kc,nr
      double precision z,a,b,enu,rofi(nr),g0(alpha),gfn(2,2,2,2),bm(nr+2),vm(nr+2),
     .  g1out(nr),g2out(nr),f1out(nr),f2out(nr),g1in(nr),g2in(nr),f1in(nr),f2in(nr),
     .  valslo(2,2,2)
C ... Local parameters
      integer n,nrk,ncouple,kr,nrmx
      double precision J,bmint,c,csq,df11,df12,df13,df14,df21,df22,df23,
     .  df24,dg11,dg12,dg13,dg14,dg21,dg22,dg23,dg24,dx,enu1,enu2,e1,e2,
     .  f1c,f2c,ff,g1c,g2c,Zbyc,mu,nuf1int,nuf2int,nug1int,nug2int,q,r,
     .  rmax,sqru,u,up,upp,vmint,kap,kapp,ug1,uf1,rg,sclg
C     For starting point
      double precision gam,gam2,X,Y,X2,Y2,gamc
C     double precision g1sin(nr),g2sin(nr),f1sin(nr),f2sin(nr)
C     double precision g1sout(nr),g2sout(nr),f1sout(nr),f2sout(nr)
C     Help arrays holding factors in the 1st order differential equation
      double precision nug1(nr),nug2(nr),nuf1(nr),nuf2(nr)
C     double precision vsnr(2,2)

      common /cc/ c

      call dpzero(g1out,nr); call dpzero(f1out,nr); call dpzero(g2out,nr); call dpzero(f2out,nr)
      call dpzero(g1in,nr); call dpzero(f1in,nr); call dpzero(g2in,nr); call dpzero(f2in,nr)
C     call dpzero(g1sout,nr); call dpzero(f1sout,nr); call dpzero(g2sout,nr); call dpzero(f2sout,nr)
C     call dpzero(g1sin,nr); call dpzero(f1sin,nr); call dpzero(g2sin,nr); call dpzero(f2sin,nr)

      csq = c*c
      Zbyc = 2*z/c
      dx = 1d0
      ff = 1d0
      rmax = rofi(nr)
      enu1 = enu; enu2 = enu ! Keep them distinct internally
      mu = imu - l - 1.5d0
      up = mu/(l-0.5d0) ; u = mu/(l+0.5d0) ; upp = mu/(l+1.5d0); sqru = dsqrt(1d0-u*u)
      e2 = enu2
      ncouple = 0
      if (mod(mode,10) == 1 .or. sqru < 1d-12) then
        sqru = 0
        ncouple = 1
        if (alpha == 2) ncouple = -1
      endif
      if (imu == 1 .or. imu == 2*l+2) then
        kap = -l-1
        up = upp
        ug1 = -upp
        uf1 = -u
        e1 = enu2
      else
        kap = l
        ug1 = +up
        uf1 = +u
        e1 = enu1
      endif
      kapp = -l-1

C ... Coefficients to the coupled differential equations; see Remarks
      do  n  = 2, nr
        nug1(n) = 1d0 + 2*z/rofi(n)/csq + ff*(e1 - vm(n) + ug1*bm(n))/csq
        nuf1(n) =    - (2*z/rofi(n)     +     e1 - vm(n) + uf1*bm(n))
        nug2(n) = 1d0 + 2*z/rofi(n)/csq + ff*(e2 - vm(n) - upp*bm(n))/csq
        nuf2(n) =    - (2*z/rofi(n)     +     e2 - vm(n) - u*bm(n))
      enddo

C --- Outwards Integration ---
      if (mod(mode/10,10) == 2) goto 10
C ... Initial points for r -> 0: V = V0 - 2Z/r, B = const
      g1out(1) = 0d0 ; f1out(1) = 0d0; g2out(1) = 0d0 ; f2out(1) = 0d0
      sclg = 1; gamc = 0
      do  n = 2, 4
        if (alpha == 1) then
          call rdinit(kap,kapp,Zbyc,bm(n),ug1,uf1,sqru,c,e1-Vm(n),gam,gam2,q,X,Y,X2,Y2)
          rg = rofi(n)**gam
          if (mod(mode/100,10) == 1) gamc = gam
          if (g0(alpha) /= 0 .and. n == 2) sclg = g0(alpha)/(rg * (1 + X*rofi(n)))
          sclg = sclg/rofi(n)**gamc
          g1out(n) = rg * (1 + X*rofi(n)) * sclg
          f1out(n) = rg * (q + Y*rofi(n)) * sclg
          g2out(n) = rg * (0 + X2*rofi(n)) * sclg
          f2out(n) = rg * (0 + Y2*rofi(n)) * sclg
        else
          call rdinit(kapp,kap,Zbyc,bm(n),ug1,uf1,sqru,c,e2-Vm(n),gam2,gam,q,X,Y,X2,Y2)
          rg = rofi(n)**gam2
          if (mod(mode/100,10) == 1) gamc = gam2
          if (g0(alpha) /= 0 .and. n == 2) sclg = g0(alpha)/(rg * (1 + X*rofi(n)))
          sclg = sclg/rofi(n)**gamc
          g1out(n) = rg * (0 + X2*rofi(n)) * sclg
          f1out(n) = rg * (0 + Y2*rofi(n)) * sclg
          g2out(n) = rg * (1 + X*rofi(n)) * sclg
          f2out(n) = rg * (q + Y*rofi(n)) * sclg
        endif
      enddo

C ... Outward integration for the first nrk points, using Runge-Kutta
      J = a*(rofi(2) + b) ; q = dexp(a/2) ; r = rofi(2)
      nrk = nr; if (kc > 0) nrk = kc  ! Integrate outward to nrk
C     print *, '!!'; nrk = nr
      n  = 2                             ! Integration begins at this point
      dg11 = 0; df11 = 0; dg21 = 0; df21 = 0
      dg12 = 0; df12 = 0; dg22 = 0; df22 = 0
      dg13 = 0; df13 = 0; dg23 = 0; df23 = 0
      g1out(1) = 0; f1out(1) = 0; g2out(1) = 0; f2out(1) = 0
      do  while (n < nrk)
        g1c = g1out(n); f1c = f1out(n); g2c = g2out(n); f2c = f2out(n)

        if (ncouple >= 0) then
          dg11 = dx*J* (c*nug1(n)*f1c                     - (kap+gamc)/r*g1c)
          df11 = dx*J* ( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + (kap-gamc)/r*f1c)
        endif
        if (ncouple <= 0) then
          dg21 = dx*J* (c*nug2(n)*f2c                     - (kapp+gamc)/r*g2c)
          df21 = dx*J* ( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + (kapp-gamc)/r*f2c)
        endif
        g1c = g1out(n) + dg11/2d0; f1c = f1out(n) + df11/2d0
        g2c = g2out(n) + dg21/2d0; f2c = f2out(n) + df21/2d0

        J = J*q; r = J/a - b

        vmint = (5*vm(n) + 2*vm(n+1) + vm(n+2))/8d0
        bmint = (5*bm(n) + 2*bm(n+1) + bm(n+2))/8d0

        nug1int = 1d0 + 2*z/r/csq + ff*(e2 - vmint + ug1 * bmint)/csq
        nuf1int =     - 2*z/r     -    (e2 - vmint + uf1 * bmint)
        nug2int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
        nuf2int =     - 2*z/r     -    (e2 - vmint - u   * bmint)

        if (ncouple > 0) then
          dg12 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          g1c = g1out(n) + dg12/2d0 ; f1c = f1out(n) + df12/2d0

          dg13 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          g1c = g1out(n) + dg13 ; f1c = f1out(n) + df13
        elseif (ncouple < 0) then
          dg22 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)
          g2c = g2out(n) + dg22/2d0 ; f2c = f2out(n) + df22/2d0

          dg23 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)
          g2c = g2out(n) + dg23 ; f2c = f2out(n) + df23
        else
          dg12 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          dg22 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)

          g1c = g1out(n) + dg12/2d0 ; f1c = f1out(n) + df12/2d0
          g2c = g2out(n) + dg22/2d0 ; f2c = f2out(n) + df22/2d0

          dg13 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          dg23 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)

          g1c = g1out(n) + dg13 ; f1c = f1out(n) + df13
          g2c = g2out(n) + dg23 ; f2c = f2out(n) + df23
        endif

        J = J*q ; r = J/a - b   ! r = rofi(n+1)

        if (ncouple >= 0) then
          dg14 = dx*J*(c*nug1(n+1)*f1c                       - (kap+gamc)/r*g1c)
          df14 = dx*J*( (nuf1(n+1)*g1c - sqru*bm(n+1)*g2c)/c + (kap-gamc)/r*f1c)
          g1out(n+1) = g1out(n) + (dg11 + 2*(dg12+dg13) + dg14)/6d0
          f1out(n+1) = f1out(n) + (df11 + 2*(df12+df13) + df14)/6d0
C         print *, 'o', n+1, f1out(n+1)/g1out(n+1), rofi(n+1) - r
        endif
        if (ncouple <= 0) then
          dg24 = dx*J*(c*nug2(n+1)*f2c                       - (kapp+gamc)/r*g2c)
          df24 = dx*J*( (nuf2(n+1)*g2c - sqru*bm(n+1)*g1c)/c + (kapp-gamc)/r*f2c)
          g2out(n+1) = g2out(n) + (dg21 + 2*(dg22+dg23) + dg24)/6d0
          f2out(n+1) = f2out(n) + (df21 + 2*(df22+df23) + df24)/6d0
        endif

        n = n + 1
      enddo

      if (gamc /= 0) then
        do  n = 2, nrk
          sclg = rofi(n)**gamc
          g1out(n) = sclg*g1out(n)
          f1out(n) = sclg*f1out(n)
          g2out(n) = sclg*g2out(n)
          f2out(n) = sclg*f2out(n)
        enddo
      endif
C      call prrmsh('g1out',rofi,g1out,nr,nr,1)
C      call prrmsh('f1out',rofi,f1out,nr,nr,1)

C ... g and derivative dg/dr at nrk
      kr = nrk; if (kc < 0) kr = -kc
      valslo(1,1,1) = g1out(kr)
      valslo(2,1,1) = g2out(kr)
      valslo(1,1,2) = (c*nug1(kr)*f1out(kr) - kap/rofi(kr)*g1out(kr))
      valslo(2,1,2) = (c*nug2(kr)*f2out(kr) - kapp/rofi(kr)*g2out(kr))

C     Equivalently
C      call rdeqvs(1,l,imu,1,1,enu,rofi(kr),z,vm(kr),bm(kr),
C     .  [g1out(kr),f1out(kr),g2out(kr),f2out(kr)],[X],vsnr)
C      print *, valslo(:,1,:) - vsnr

C --- Inwards Integration ---
   10 continue
      if (kc == nr .or. mod(mode/10,10) == 1) goto 20

      r = rofi(nr)
      g1in(nr) = gfn(1,1,1,alpha)
      f1in(nr) = gfn(2,1,1,alpha)
      g2in(nr) = gfn(1,2,1,alpha)
      f2in(nr) = gfn(2,2,1,alpha)
      gamc = 0

C ... Runge-Kutta for points nr to nrk
      nrmx = nr

      nrk = 2; if (kc > 0) nrk = kc ! Integrate inward nrmx to nrk

C     For debugging
C      print *, '!!'
C      nrmx = 777; nrk = 2
C      g1in(nrmx) = g1out(nrmx)
C      g2in(nrmx) = g2out(nrmx)
C      f1in(nrmx) = f1out(nrmx)
C      f2in(nrmx) = f2out(nrmx)


      J = a*(rofi(nrmx) + b) ; q = dexp(-a/2d0) ; r = rofi(nrmx)

      do  n = nrmx, nrk+1, -1
        g1c = g1in(n); f1c = f1in(n); g2c = g2in(n); f2c = f2in(n)

        if (ncouple >= 0) then
          dg11 = dx*J* (c*nug1(n)*f1c                     - (kap+gamc)/r*g1c)
          df11 = dx*J* ( (nuf1(n)*g1c - sqru*bm(n)*g2c)/c + (kap-gamc)/r*f1c)
        endif
        if (ncouple <= 0) then
          dg21 = dx*J* (c*nug2(n)*f2c                     - (kapp+gamc)/r*g2c)
          df21 = dx*J* ( (nuf2(n)*g2c - sqru*bm(n)*g1c)/c + (kapp-gamc)/r*f2c)
        endif
        g1c = g1in(n) - dg11/2d0; f1c = f1in(n) - df11/2d0
        g2c = g2in(n) - dg21/2d0; f2c = f2in(n) - df21/2d0

        J = J*q; r = J/a - b

        vmint = (5*vm(n) + 2*vm(n-1) + vm(n-2))/8d0
        bmint = (5*bm(n) + 2*bm(n-1) + bm(n-2))/8d0

        nug1int = 1d0 + 2*z/r/csq + ff*(e2 - vmint + ug1 * bmint)/csq
        nuf1int =     - 2*z/r     -    (e2 - vmint + uf1 * bmint)
        nug2int = 1d0 + 2*z/r/csq + ff*(e2 - vmint - upp * bmint)/csq
        nuf2int =     - 2*z/r     -    (e2 - vmint - u   * bmint)

        if (ncouple > 0) then
          dg12 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          g1c = g1in(n) - dg12/2d0 ; f1c = f1in(n) - df12/2d0

          dg13 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          g1c = g1in(n) - dg13 ; f1c = f1in(n) - df13
        elseif (ncouple < 0) then
          dg22 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)
          g2c = g2in(n) - dg22/2d0 ; f2c = f2in(n) - df22/2d0

          dg23 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)
          g2c = g2in(n) - dg23 ; f2c = f2in(n) - df23
        else
          dg12 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df12 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          dg22 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df22 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)

          g1c = g1in(n) - dg12/2d0 ; f1c = f1in(n) - df12/2d0
          g2c = g2in(n) - dg22/2d0 ; f2c = f2in(n) - df22/2d0

          dg13 = dx*J*(c*nug1int*f1c                     - (kap+gamc)/r*g1c)
          df13 = dx*J*( (nuf1int*g1c - sqru*bmint*g2c)/c + (kap-gamc)/r*f1c)
          dg23 = dx*J*(c*nug2int*f2c                     - (kapp+gamc)/r*g2c)
          df23 = dx*J*( (nuf2int*g2c - sqru*bmint*g1c)/c + (kapp-gamc)/r*f2c)

          g1c = g1in(n) - dg13 ; f1c = f1in(n) - df13
          g2c = g2in(n) - dg23 ; f2c = f2in(n) - df23
        endif

        J = J*q ; r = J/a - b   ! r = rofi(n-1)

        if (ncouple >= 0) then
          dg14 = dx*J*(c*nug1(n-1)*f1c                     - (kap+gamc)/r*g1c)
          df14 = dx*J*( (nuf1(n-1)*g1c - sqru*bm(n-1)*g2c)/c + (kap-gamc)/r*f1c)
          g1in(n-1) = g1in(n) - (dg11 + 2*(dg12+dg13) + dg14)/6d0
          f1in(n-1) = f1in(n) - (df11 + 2*(df12+df13) + df14)/6d0
C         print *, 'i', n-1, f1in(n-1)/g1in(n-1), rofi(n-1) - r
        endif
        if (ncouple <= 0) then
          dg24 = dx*J*(c*nug2(n-1)*f2c                     - (kapp+gamc)/r*g2c)
          df24 = dx*J*( (nuf2(n-1)*g2c - sqru*bm(n-1)*g1c)/c + (kapp-gamc)/r*f2c)
          g2in(n-1) = g2in(n) - (dg21 + 2*(dg22+dg23) + dg24)/6d0
          f2in(n-1) = f2in(n) - (df21 + 2*(df22+df23) + df24)/6d0
        endif

      enddo

C      call prrmsh('g1in',rofi,g1in,nrmx,nrmx,1)
C      call prrmsh('f1in',rofi,f1in,nrmx,nrmx,1)

C ... g and derivative dg/dr at nrk
      kr = nrk; if (kc < 0) kr = -kc
      valslo(1,2,1) = g1in(kr)
      valslo(2,2,1) = g2in(kr)
      valslo(1,2,2) = (c*nug1(kr)*f1in(kr) - kap/rofi(kr)*g1in(kr))
      valslo(2,2,2) = (c*nug2(kr)*f2in(kr) - kapp/rofi(kr)*g2in(kr))

C     Equivalently
C      call rdeqvs(1,l,imu,1,1,enu,rofi(kr),z,vm(kr),bm(kr),
C     .  [g1in(kr),f1in(kr),g2in(kr),f2in(kr)],[X],vsnr)
C      print *, valslo(:,2,:) - vsnr

   20 continue

C      call info5(10,1,-1,' rdeqinout exit l=%i, imu=%i, mu=%d, alpha=%i',l,imu,mu,alpha,0)
C      call info5(10,0,-1,' out g1 f1 g2 f2(nr) = %g %g %g %g',g1out(nr),f1out(nr),g2out(nr),f2out(nr),0)
C      call info5(10,0,0,'  in g1 f1 g2 f2(nr) = %g %g %g %g',g1in(nr),f1in(nr),g2in(nr),f2in(nr),0)

      end subroutine rdeqinout

      subroutine gftopv(cpin,imu,l,alpha,ampf,ampg,gfn,Pv)
C- Copy and scale components of Dirac function at one point to/from P vector
C ----------------------------------------------------------------------
Ci Inputs
Ci   cpin  : T, copy/scale gf to P;  F, copy/scale P to gf
Ci   imu   :
Ci   l     :
Ci   alpha : 1 or 2
Ci   ampf  :
Ci   ampg  :
Ci   gfn   :
Ci   Pv    :
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   22 Oct 14
C ----------------------------------------------------------------------
      implicit none
      logical cpin
      integer imu,l,alpha
      double precision ampf,ampg,gfn(2,2,2,2),pV(4)

      if ((imu == 1 .or. imu == 2*l+2) .and. alpha == 1) then
        if (cpin) then
          Pv(1) = gfn(2,1,1,alpha)/gfn(1,1,1,alpha)*ampf ! value of f1 at the MT boundary
        else
!         gfn(2,1,1,alpha) = Pv(1)/ampf
          gfn(2,1,1,alpha) = Pv(1)*gfn(1,1,1,alpha)/ampf
        endif
      elseif ((imu == 1 .or. imu == 2*l+2) .and. alpha == 2) then
        if (cpin) then
          Pv(1) = gfn(2,2,1,alpha)/gfn(1,2,1,alpha)*ampf ! value of f2 at the MT boundary
        else
!         gfn(2,2,1,alpha) = Pv(1)/ampf
          gfn(2,2,1,alpha) = Pv(1)*gfn(1,2,1,alpha)/ampf
        endif
C      elseif ((f1out2(nr)==0).and.(f2out(nr)==0)) then
C        if (alpha == 1) then
C          gfn(1,1,1,alpha) = Pv(1)
C          gfn(2,1,1,alpha) = Pv(2)
C        else
C         gfn(1,2,1,alpha) = Pv(3)
C         gfn(2,2,1,alpha) = Pv(4)
C       endif
      else
        if (cpin) then
          Pv(1) =  gfn(1,1,1,alpha)*ampg
          Pv(2) =  gfn(2,1,1,alpha)*ampf
          Pv(3) =  gfn(1,2,1,alpha)*ampg
          Pv(4) =  gfn(2,2,1,alpha)*ampf
        else
          gfn(1,1,1,alpha) = Pv(1)/ampg
          gfn(2,1,1,alpha) = Pv(2)/ampf
          gfn(1,2,1,alpha) = Pv(3)/ampg
          gfn(2,2,1,alpha) = Pv(4)/ampf
        endif
      endif
      end

      subroutine prodrw(nr,ncomp,g1,g2,wt,normn)
C- Vector inner product of n-component function, given quadrature weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of mesh points
Ci   ncomp :number of components in each vector
Ci   g1    :First  n-vector
Ci   g2    :Second n-vector
Co   wt    :weights for numerical integration (radwgt)
Co Outputs
Ci   normn :result of inner product
Cr Remarks
Cr
Cu Updates
Cu   19 Mar 15 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,ncomp
      double precision normn,g1(nr,ncomp),g2(nr,ncomp),wt(nr)
C ... Local parameters
      integer i
      procedure(real(8)) :: dot3

      normn = 0
      do  i = 1, ncomp
        normn = normn + dot3(nr,wt,g1(1,i),g2(1,i))
      enddo
      end subroutine prodrw

      real(8) function diracp(nr,w1,w2,wt)
C- Inner product of a pair of Dirac 4-vectors on a radial mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   w1    :First 4-vector
Ci   w2    :Second 4-vector
Co   wt    :weights for numerical integration (radwgt)
Ci   a     :Parameter defining radial mesh; see radmsh.f
Ci   b     :radial mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Co Outputs
Ci   normn :result of inner product
Cu Updates
Cu   29 Sep 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr
      double precision w1(nr,4),w2(nr,4),wt(nr)
C ... Local parameters
      integer i
      procedure(real(8)) :: dot3

      diracp = 0
      do  i = 1, 4
        diracp = diracp + dot3(nr,wt,w1(1,i),w2(1,i))
      enddo
C     print *, 'norm',diracp, dot3(4*nr,wt,w1,w2)

      end

      subroutine maxsearch(g,nr,maxp)
C- Find the first max going from nr backwards
      implicit none
C ... Passed parameters
      integer nr
      double precision g(nr)
      integer, intent(out) :: maxp
C ... Local parameters
      double precision D0,D1,sig
      integer n

      D0 = g(nr) - g(nr-1)
      if (D0 == 0) then
        maxp = 1
        return
      endif
      do  n = nr-1, 2, -1
        D1 = g(n) - g(n-1)
        maxp = n
        if ((D1 == 0).or.(D0 == 0)) return
        sig = D1/D0
        if (sig <= 0) return
        D0 = D1
      enddo

      end subroutine maxsearch

      subroutine rdeqvs(mode,l,imu,kr,nr,enu,rofi,z,vm,bm,psi,gfn,valslo)
C- Return value and slope of large components of Dirac vector at point kr
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit
Ci         :0 make valslo from gfn
Ci         :1 like mode=0, but take val,slo from psi(kr,:)
Ci         :2 return valslo, gfn from data stored in prior call
Ci         :10s digit
Ci         :0 forward copy
Ci         :1 reverse the sense of the copy (psi or gfn is returned)
Ci         :  This digit isn't relevant if 1s digit mode = 2
Ci   l     :l quantum number
Ci   imu   :projection of j on z axis (positive integer, counts -l-1/2...l+1/2)
Ci   kr    :Point at which to return value and slopes
Ci   nr    :number of radial mesh points
Ci   enu   :enu's for making charge density
Ci   rofi  :radial mesh points
Ci   z     :nuclear charge
Ci   vm    :scalar potential
Ci   bm    :magnetic field
Cio Inputs/Outputs
Cio  psi   :psi(nr,:) contains same information as gfn
Cio  gfn   :Function values at rofi(kr)
Cio        :gfn(1) = g1(kr)   gfn(2) = f1(kr)    gfn(3) = g2(kr)   gfn(4) = f2(kr)
Cio  valslo: values and slopes of large components gi at matching point kc
Cio        : valslo(i,k) : i=1,2 labels kappa; k=1,2 for (val,slo)
Cio        : valslo(i,1)   valslo(i,2)
Cio        :    g_i           g'_i
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   21 Jul 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,l,imu,kr,nr
      double precision enu,z,rofi(nr),psi(nr,4),vm(nr),bm(nr),gfn(4),valslo(2,2)
C ... Local parameters
      integer mode0
      double precision c,csq,nug1,nug2,kap,kapp,mu,u,up,upp,gfnl(4)
      real(8),parameter:: NULLR =-99999
      real(8) :: savevs(2,2) = NULLR
      real(8) :: savegf(4) = NULLR
      save savevs, savegf
      common /cc/ c   ! Speed of light, or infinity in nonrelativistic case

      mode0 = mod(mode,10)
      gfnl = gfn

      kap = l; kapp = -l-1; mu  = imu - l - 1.5d0; csq = c*c
      if (imu == 1 .or. imu == 2*l+2) kap = -l-1
      up = mu/(l-0.5d0) ; u = mu/(l+0.5d0) ; upp = mu/(l+1.5d0)
      if (imu == 1 .or. imu == 2*l+2) then
        nug1 = 1d0 + 2*z/rofi(kr)/csq + (enu - vm(kr) - upp*bm(kr))/csq
      else
        nug1 = 1d0 + 2*z/rofi(kr)/csq + (enu - vm(kr) + up *bm(kr))/csq
      endif
      nug2 = 1d0 + 2*z/rofi(kr)/csq + (enu - vm(kr) - upp*bm(kr))/csq

C     If from prior call, copy savevs to valslo and savegf to gfnl
      if (mode0 == 2) then
        call dcopy(4,savevs,1,valslo,1)
        call dcopy(4,savegf,1,gfn,1)
        return
      elseif (mode0 == 1) then
        gfnl(:) = psi(kr,:)
        mode0 = 0
      endif

C     Forward copy : make valslo from gfnl
      if (mod(mode/10,10) == 0) then
        valslo(1,1) = gfnl(1)
        valslo(2,1) = gfnl(3)
        valslo(1,2) = c*nug1*gfnl(2) - kap/rofi(kr)*gfnl(1)
        valslo(2,2) = c*nug2*gfnl(4) - kapp/rofi(kr)*gfnl(3)
C     Reverse copy : make gfnl from valslo; also copy to psi(kr,:)
      else
        gfnl(1) = valslo(1,1)
        gfnl(3) = valslo(2,1)
        gfnl(2) = (rofi(kr)*valslo(1,2) +  kap*valslo(1,1))/(c*nug1*rofi(kr))
        gfnl(4) = (rofi(kr)*valslo(2,2) + kapp*valslo(2,1))/(c*nug2*rofi(kr))
        if (mode0 == 1) psi(kr,:) = gfnl(:)
        if (mode0 == 0) gfn = gfnl
      endif
      savegf = gfnl
      savevs = valslo

      end

      subroutine rdinit(kap1,kap2,Zbyc,Bf,ug1,uf1,sqru,c,EmV,
     .  gam1,gam2,q,X,Y,X2,Y2)
C- Parameters for Dirac equation as r->0
C ----------------------------------------------------------------------
Ci Inputs
Ci   kap1  :interstitial kinetic energy
Ci   kap2  :interstitial kinetic energy
Ci   kap2  :kinetic energies for which strux are calculated
Ci   Zbyc  :2*Z/c
Ci   Bf    :Magnetic field for small r
Ci   ug1   :Coefficient to magnetic field coupling term; see Remarks
Ci   uf1   :Coefficient to magnetic field coupling term; see Remarks
Ci   sqru  :Coefficient to magnetic field coupling term; see Remarks
Ci   c     :Speed of light
Co Outputs
Co   gam1  :sqrt(kap1**2 - Zbyc**2)
Co   gam2  :sqrt(kap2**2 - Zbyc**2)
Co   q     :(gam1+kap1)/Zbyc
Co   EmV   :E - V(0), where E is the energy
Co   X     :Coefficient to term r^(gam1+1) in g1
Co   Y     :Coefficient to term r^(gam1+1) in f1
Co   X2    :Coefficient to term r^(gam1+1) in g2
Co   Y2    :Coefficient to term r^(gam1+1) in f2
Cl Local variables
Cl   B1    :
Cl   B2    :
Cl   B4    :
Cr Remarks
Cr   Dirac equation for r -> 0: V = V(0) - 2Z/r, B = const
Cr   Starting from the Dirac equations at small r written as:
Cr   (d/dr + kap1/r)*g1(r) - ((E-V(0)+(c*Zbyc)/r)/c^2 + 1 + B1) * c*f1(r)
Cr c*(d/dr - kap1/r)*f1(r) +  (E-V(0)+(c*Zbyc)/r + B2) * g1(r)
Cr   (d/dr + kap2/r)g2(r)  - (c*Zbyc)/r/c^2                     * c*f2(r)
Cr c*(d/dr - kap2/r)f2(r)  + (c*Zc)/r * g2(r) + B4*g1(r)
Cr   g can be expanded as a polynomial
Cr     g1 = r^gam1*(1 + X*r)
Cr     f1 = r^gam1*(q + Y*r)
Cr     g2 = r^gam1*(0 + X2*r)
Cr     f2 = r^gam1*(0 + Y2*r)
Cr   where
Cr     gam1^2 = kap1**2 - Zbyc**2
Cr          q = (gam1+kap1)/Zbyc
Cr       Zbyc = 2*Z/c
Cr     X and Y and X2 and Y2 are long expressions, given below
Cr     Turek's 6.19 defines a, which is gam1-1.  Note  g(rdeq) = r*g(Turek)
Cu Updates
Cu   15 Jul 16
C ----------------------------------------------------------------------
      implicit none

      double precision kap1,kap2,Zbyc,ug1,uf1,sqru,c,
     .  gam1,gam2,q,EmV,B1,B2,B4,X,Y,X2,Y2
      double precision Bf(*)

      gam1 = dsqrt(kap1**2 - Zbyc**2)
      gam2 = dsqrt(kap2**2 - Zbyc**2)
      q = (gam1+kap1)/Zbyc
      B1 = ug1*Bf(2)/c**2 ! Magnetic factor in coefficient to +f1, Eq. (D1)
      B2 = uf1*Bf(2)      ! Magnetic factor in coefficient to -g1, Eq. (D2)
      B4 = sqru*Bf(2)/c   ! Magnetic factor in coefficient to -g1, Eq. (D4)
      X = (c**2*Zbyc**2 + (kap1+gam1-2*Zbyc**2)*(EmV+c**2)
     .  -  B1*c**2*(Zbyc**2-kap1-gam1) - B2*Zbyc**2) / (c*(2*gam1+1)*Zbyc)
      Y = ((-2*kap1*EmV)-2*gam1*EmV-EmV-c**2*kap1-c**2*gam1
     .    - B1*c**2*(kap1+gam1) - B2*(kap1+gam1+1)) / (c*(2*gam1+1))
      X2 = (Zbyc*B4)/(c*(gam2-gam1-1)*(gam2+gam1+1))
      Y2 = (kap2+gam1+1)/Zbyc*X2

      end

      subroutine orthdirac(opt,nr,wt,psi1,psi2)
C- Orthogonalize Dirac psi1 to psi2
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 Do not orthogonalize
Ci         :1 Orthogonalize psi2 to psi1
Ci         :  psi1 is unchanged unless it is renormalized.
Ci         :  Recommended when psi1 and psi2 are solutions of different equations
Ci         :  esp psi1 is w.f and psi2 is energy derivative
Ci         :2 Orthogonalize psi1 and psi2 in a symmetric manner; See Remarks.
Ci         :  psi1 and psi2 are both modified; the norm of neither is conserved.
Ci         :  Recommended when psi1 and psi2 are two independent solutions
Ci         :  of the same equation.
Ci         :10's digit
Ci         :1  normalize psi1 and psi2  after possible orthogonalization
Ci   nr    :number of radial mesh points
Ci   wt    :radial mesh weights
Cio Inputs/Outputs
Ci   psi1  :first Dirac 4-vector
Ci   psi2  :Second Dirac 4-vector
Co Outputs
Co   s
Cr Remarks
Cr   psi1 and psi2 are a pair of independent solutions to the Dirac equation
Cr   or an energy derivative.
Cr   This routine orthogonalizes psi1 and psi2 orthogonalized to each other
Cr   Prescription for symmetric orthogonalization (opt=1)
Cr     Let P1 and P2 be normalized psi1 and psi2.
Cr     Define
Cr       P1 = psi1/sqrt(s11); P2 = psi2/sqrt(s22)
Cr       P1' = P1 - xx M12 P2
Cr       P2' = P2 - xx M12 P1
Cr       <P1',P2'> = 0
Cr     provided
Cr       xx = -(sqrt(1-M12^2)-1)/M12^2
Cr       M12 = (P1,P2) = 1/sqrt(s11 s22) s12
Cr       sij = (psii,psij)
Cr     Then
Cr       psi1' = sqrt(s11) P1' = psi1 - sqrt(s11) xx M12 P2
Cr                             = psi1 - xx s12/s22 psi2
Cr       psi2' = sqrt(s22) P2' = psi2 - sqrt(s22) xx M12 P1
Cr                             = psi2 - xx s12/s11 psi1
Cr       <psi1'|psi2'> = 0
Cu Updates
Cu   02 Aug 16
C ----------------------------------------------------------------------
      implicit none
      integer opt,nr
      double precision wt(nr),psi1(nr,4),psi2(nr,4)
C ... Local parameters
      double precision xx,m12
      double precision s11,s12,s22,p1(nr,2,2),p2(nr,2,2)
      real(8), parameter :: tol = 1d-12
      procedure(real(8)) :: dot3,diracp,dlength

      s12 = diracp(nr,psi1,psi2,wt)
      if (abs(s12) < tol) return
      s11 = diracp(nr,psi1,psi1,wt)
      s22 = diracp(nr,psi2,psi2,wt)

C     debugging
C      write(*,"(/'% rows 2 cols 2  overlap on input')")
C      print 333, s11, s12
C      print 333, s12, s22
C  333 format(2f18.12)

C --- Orthogonalize ---
      if (mod(opt,10) == 1) then
        call daxpy(nr*4,-s12/s11,psi1,1,psi2,1)
      elseif (mod(opt,10) == 2) then
        call dpcopy(psi1,p1,1,4*nr,1d0)
        call dpcopy(psi2,p2,1,4*nr,1d0)
        m12 =  s12/dsqrt(s11*s22)
        xx = 0.5d0 + M12**2/8
        if (abs(M12) > 2d-3) then
          xx = -(dsqrt(1-M12**2)-1)/M12**2
        endif
        call daxpy(nr*4,-xx*s12/s22,p2,1,psi1,1)
        call daxpy(nr*4,-xx*s12/s11,p1,1,psi2,1)
      endif

C --- Reormalize ---
      if (mod(opt/10,10) /= 0) then
        xx = diracp(nr,psi1,psi1,wt)
        call dscal(4*nr,1/dsqrt(xx),psi1,1)
        xx = diracp(nr,psi2,psi2,wt)
        call dscal(4*nr,1/dsqrt(xx),psi2,1)
      endif

C     debugging
C      s11 = diracp(nr,psi1,psi1,wt)
C      s12 = diracp(nr,psi1,psi2,wt)
C      s22 = diracp(nr,psi2,psi2,wt)
C      write(*,"(/'% rows 2 cols 2  overlap on output')")
C      print 333, s11, s12
C      print 333, s12, s22

      end

C      subroutine snot(i,debug)
C      logical debug
C      integer i
CC     return
C      print *, 'snot',i
C      call setpr(100)
C      debug = .true.
C      end
C      subroutine snit(iter,debug,gfn)
C      integer iter
C      logical debug
C      double precision gfn(2,2,2,2)
C      character*100 outs
C      integer a2vec,i,ix(10)
C
C      if (.not. debug) return
C
CC      return
C
C      if (iter >= 14 .or. iter < 0 .or. .true.) then
C        print *, 'gfn=',gfn(:,1,1,1),'  change to?'
C        read(*,'(a100)') outs
C        if (outs == ' ') return
C        i = 0
C        i = a2vec(outs,len(outs),i,4,', ',2,3,2,ix,gfn(:,1,1,1))
C        iter = -1
C      endif
C      end
C      subroutine snut(kc,nre,psi)
C      integer kc,nre
C      double precision psi(nre,4)
C      print *, psi(kc-1,:)
C      print *, psi(kc,:)
C      print *, psi(kc+1,:)
C      end

      subroutine vpxtrap(nr,nxtrap,npoly,a,b,rofi,v,errmx)
C- Extrapolate potential beyond the MT radius
      implicit none
      integer nr,nxtrap,npoly
      double precision a,b,errmx,rofi(nr),v(nr+nxtrap)
      integer i
      double precision rmax

C     Extrapolate V(r) by 2 radial mesh points to allow Runge Kutta integration to nr
      do  i = 1, nxtrap
        rmax = exp(i*a)*(rofi(nr)+b)-b   ! Equivalently:  rmax = b*(dexp((nr-1+i)*a)-1) = rofi(nr+i)
        call polinx(rofi(nr-4),v(nr-npoly+1),npoly,rmax,1d-6,v(nr+i),errmx)
      enddo
      end

      subroutine testdirac(nr,l,imu,E,Z,rofi,psi,v)
C- Check precision for a solution of a 2-component Dirac equation (B=0)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   l
Ci   imu
Ci   E
Ci   z     :nuclear charge
Ci   rofi  :radial mesh points
Ci   psi   :r*g and r*f, g and f = large and small components of Dirac eqn
Ci   v     :spherical potential (atomsr.f)
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  From Turek's book, Eq. 6.13
Cr    (d/dr + (1+k)/r) gtrue_k - (1+(E-V(r))/c^2) c*ftrue_k = 0
Cr    (d/dr + (1-k)/r) ftrue_k +   ((E-V(r)))/c     gtrue_k = 0
Cr  Note psi(r) = r(g,f)true = (g,f)
Cr  Use : r * d/dr (gtrue_k) = d/dr g_k - g_k/r
Cr  Substituting (g,f) for (gtrue,ftrue):
Cr    (d/dr + k/r) g_k - (1+(E-V(r))/c^2) c*f_k = 0
Cr    (d/dr - k/r) f_k +    (E-V(r))/c g_k      = 0
Cu Updates
Cu   18 Aug 16
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,l,imu
      double precision E,Z,rofi(nr,2),psi(nr,2),v(nr)
C ... Local parameters
      logical, parameter :: lrat=.false.
      integer ic,lerr,kr,kap
      integer :: npoly=12
      double precision psip(nr,2)
      integer, parameter :: NULLI=-99999
      real(8) :: mu,nug,c,csq,err1,err2,errint1,errint2,tol=1d-10

      common /cc/ c

C     call prrmsh('psi',rofi,psi,nr,nr,2)

C     Compute dg/dr and df/dr
      do  ic = 1, 2
        call poldvm(rofi,psi(1,ic),nr,npoly,lrat,tol,lerr,psip(1,ic))
      enddo

C     call prrmsh("psi'",rofi,psip,nr,nr,2)

      mu = imu - l - 1.5d0; csq = c*c
      if (imu == 1 .or. imu == 2*l+2) then
        kap = -l-1
C        up = upp
C        ug1 = -upp
C        uf1 = -u
C        e1 = enu2
      else
        kap = l
C        ug1 = +up
C        uf1 = +u
C        e1 = enu1
      endif

   10 continue
      errint1 = 0; errint2 = 0
      do  kr = 2, nr

        nug = 1d0 + (E + 2*Z/rofi(kr,1) - V(kr))/csq
        err1 = psip(kr,1) + kap/rofi(kr,1)*psi(kr,1) - nug*c*psi(kr,2)
        err2 = (psip(kr,2) - kap/rofi(kr,1)*psi(kr,2))*c + (E + 2*Z/rofi(kr,1) - V(kr))*psi(kr,1)
        errint1 = errint1 + err1**2*rofi(kr,2)
        errint2 = errint2 + err2**2*rofi(kr,2)
        print *, kr, err1, err2
      enddo

      call info5(1,0,0,' E=%;12F  RMS error = %;12F %;12F ',E,errint1,errint2,4,5)

      E = -NULLI
      print *, 'E=?'; read(*,*) E
      if (E == NULLI) stop
      goto 10

      end
