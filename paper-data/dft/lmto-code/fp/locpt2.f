      subroutine locpt2(job,z,rmt,rg,a,nr,nsp,nspc,cofg,cofh,ceh,rfoc,lfoc,
     .  nlml,qmom,vval,rofi,rwgt,rho1,rho2,rhoc,nbf,bscal,bloc,rhobg,vext,
     .  rhol1,rhol2,v1,v2,v1es,v2es,s_atparms,gpotb,efg,buf)
C- Makes the potential at one site, and associated energy terms.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_atparms
Ci     Elts read:  qv1 qv2 sgpotb sgpote vales1 vales2 rmu1 rmu2 rep1
Ci                 rep2 focexc focvxc rvvxc rveps rhovv1 rhovv2 rhvxt1
Ci                 rhvxt2 rhoexc rhovxc qloc alocc atrue aloc qlocc
Co     Stored:     qv1 qv2 qlocc alocc qloc qcor1 aloc bfam atrue
Co                 sgpotb sgpote rvext vales1 vales2 vvesat cpnves
Co                 rvvxcv rvepsv rvexv rvecv rvvxc rveps rep1 rmu1 rep2
Co                 rmu2 rhovv1 rhovv2 rhvxt1 rhvxt2 rhoexc rhoex rhoec
Co                 rhovxc rhcvef1
Co     Allocated:  *
Cio    Elts passed: aloc sgpote rep1 rmu1 rep2 rmu2 focexc focex focec
Cio                focvxc rhoexc rhovxc vales1 vales2 rhovv1 rhovv2 qv1
Cio                qv2 atrue
Cio    Remarks
Cio
Ci Inputs
Ci   job   :1s to 100s digit not used here.
Ci         :1000s digit
Ci         :1 Generate valence-only rvepsv,rvexv,rvecv,rvvxcv,rveps,rvvxc
Ci         :100000s digit
Ci         : 0 use density as is.
Ci         : 1 reset points of negative density to positive density
Ci             for purposes of calculating vxc
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci         : 3 like 1, but modifed rho is retained.
Ci         :1000000s digit
Ci         : 0 normal mode
Ci         : 1 Add external potential vext to Hartree+vxc
Ci   z     :nuclear charge
Ci   rmt   :augmentation radius
Ci   rg    :smoothing radius for compensating gaussians used to
Ci         :correct the multipole moments of local smooth density
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci         :cofg = Y0 * qcorg
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   lfoc  :switch specifying treatment of core density.
Ci         :0 => core w.f. have val,slo = 0 at rmt
Ci         :1 => core included explicitly with valence
Ci         :2 => core included perturbatively
Ci   nlml  :L-cutoff for rho1,rho2
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   vval  :boundary condition for estat potential at MT boundary.
Ci         :See remarks.
Ci   rofi  :radial mesh points for density and potential
Ci   rwgt  :radial mesh weights to integrate function on radial mesh
Ci         :Integral f(r) dr = sum_i f_i * wt_i
Ci   rho1  :local valence density = sum_ilm rho1(ilm) * r**2  Y_L(ilm),
Ci         :on radial mesh.
Ci   rho2  :local smoothed valence density, defined as rho1
Ci         :Local atomic valence density is rho1-rho2
Ci   rhoc  :core density times 4*pi*r*r
Ci   nbf   :dimensions bloc
Ci   bscal :scale vxc+ - vxc- from true density by bscal(1)
Ci         :scale vxc+ - vxc- from smooth density by bscal(2)
Ci         :Relevant for nsp=2 only.  bscal(1:2) should normally be 1.
Ci   bloc  :external magnetic field for this site
Ci   rhobg:compensating background density
Ci   vext :external potential, 1-center expansion
Ci   buf   :Buffer to contain redirected standard output.
Ci         :This routine uses awrite to write output to stdout.
Ci         :It gets redirected to buf if buffer mode is turned on
Co Outputs
Co   rhol1 :full true electron density, rho1 + rhoc  (density * r**2)
Co         :Note: spin part of rhol1 is not rotated to z axis now
Co   rhol2 :full smooth density, i.e. uncompensated rho2 plus
Co         :compensating gaussians + pseudocore charge
Co         :rho2 + gval + gcor + gnuc
Co         :Note: spin part of rhol2 is not rotated to z axis now
Co   v1    :true potential at this site, excluding nuclear contribution
Co         :but including background contribution; see Remarks
Co         :See Remarks concerning boundary condition at rmt
Co         :... If vext is present, it added to v1
Co   v2    :Total smooth potential including background potential
Co         : = ves[rhol2] + ...
Co         :   ... lfoc=0 : vxc(rho2)
Co         :       lfoc=1 : vxc(rho2+sm-rhoc)
Co         :       lfoc=2 : vxc(rho2) + dvxc/drho * sm-rhoc
Co         :Apart from differences in l-truncation in the XC potential
Co         :(i.e. true and local sm. densities are different),
Co         :v1(nr,1,1)-2*z/y0/rmt = v2(nr,1,1)
Co         :See Remarks concerning boundary condition at rmt
Co         :NB: only correct so far for collinear, l-independent bfield
Co         :... If vext is present, it added to v2
Co   s_atprms : structure with the following:
Co
Co   vvesat:integral of valence density times electrostatic potential
Co         : = int rho1 * Ves1 - rho2 * Ves2
Co         :   where rho1 is valence density while Ves1 includes core rho and 2*Z/r
Co         :   and   rho2 sm valence density while Ves2 = ves[rho2+gval+gcor+gnuc]
Co          Total electrostatic energy has a term vvesat/2
Co   cpnves:integral of core+nucleus times electrostatic potential
Co          = int [ (rhoc-z) ves1~ - (rhocsm-rhonsm) ves2~ ]
Co          where ves1~ is the estat potential of the true density
Co          where ves2~ is the estat potential of the smooth density
Co          Total electrostatic energy has a term cpnves/2
Co   rhoexc:integral of valence density times xc energy density
Co   rhoex :integral of valence density times exch. energy density
Co   rhoec :integral of valence density times corr. energy density
Co   rhovxc:integral of valence density times xc potential
Co   rvepsv:integral of valence density exc(valence density)
Co   rvexv :integral of valence density ex(valence density)
Co   rvecv :integral of valence density ec(valence density)
Co   rvvxcv:integral of valence density vxc(valence density)
Co   rveps :
Co   rvvxc :
Co   rhovv1:integral valence density times full potential
Co         : = int n1_val * v1~   where v1~ = ves[n1] + vxc[n1]
Co         :used to compute double-counting for kinetic energy.
Co   rhovv2:integral valence sm-density times pseudo-potential:
Cl         : = int (n2*v2~) +  sum_L qmom_L gpotb_L
Cl         :   [- perturbation (focvxc(1)+focvxc(2)) if lfoc=2]
Co   rhcv1 :integral rhoc*(v1-2Z/r)
Co   qloc  :total valence charge in sphere rho1 - rho2
Co   qlocc :total core charge in sphere qcor1 - qcor2
Co   atrue :net magnetic moment of true density (vector)
Co   aloc  :total valence magnetic moment in sphere
Co   alocc :total core magnetic moment in sphere
Co   gpotb :gpotb(ilm) = integral [gaussians(r,radius rg)_ilm * smooth ves
Co         :gpotb = local analog of gpot0 generated by smves.
Co         :smves = v2 = ves(rhol2), w/ b.c. v2(rmax)=vval
Co         :If lvext=T, vext is added to smves.
Co         :vext = external potential, represented in a Ylm expansion.
Co   focexc:(lfoc=2): integral rhochs * vxc
Co         :otherwise, 0
Co   focex :integral of smoothed core and exchange energy
Co   focec :integral of smoothed core and correlation energy
Co   focvxc:If lfoc=2, integral rhochs * (dvxc/drho2 * rho2)
Co         :Otherwise, 0
Co   bfam  :sum (bext * local moment)/2 --- for double counting.
Co   efg   :l=2 potential at nucleus
Cl Local variables
Cl   rhochs:Part of smooth core density contains Hankel tails
Cl         :Used to make vxc[n2] and exc[n2]
Cl         :rhochs = srfpi*cofh*xi(0)*r*r
Cl   qcor1 :true core charge inside rmt
Cl   qcor2 :smoothed core charge inside rmt
Cl   rhonsm:nuclear density smoothed into gaussian of width rg
Cl   rhocsm:core density smoothed into gaussian of width rg
Cl   rvext :rvext(1,:) true rho*vext; rvext(2,:) sm rho*vext
Cl         :rvext(:,1) valence*vext; rvext(:,2) core*vext; rvext(:,3) nuc*vext
Cl         :Note: contribution from compensating gaussians could be be REMOVED from
Cl         :rvext(2,1) and the corresponding interstitial contribution (smvxt)
Cl         :For now they are retained.
Cr Remarks
Cr   Electrostatic potential from density rho confined to sphere
Cr   of radius R. Here rho(r) = rho_L(r) = n_L(r) * r**2
Cr     vL(ir) = 8*pi*(f1(ir) + f2(ir))
Cr         f1 = J(r) int_0^r [H(R)-H(r')] rho(r')
Cr         f2 = H(r) int_r^R J(r') rho(r')
Cr   Valid for any L and TF screening energy E; see routine lscoul.
Cr   H and J are Hankel, Bessel functions, MSM definitions (besslr)
Cr   No screening here (E=0), so H=Y_0/r and J=Y_0 for the l=0 case.
Cr   Boundary conditions:  at r=R, f1=0 => v(R) = 2*Q/R
Cr
Cr   On boundary conditions for the estat potential vval at the MT
Cr   boundary.  In the original formulation (see Springer book) the
Cr   the total hamiltonian and energy did not depend on the choice of
Cr   vval, as the contribution from true and smooth potentials exactly
Cr   cancel.  However, this is no longer true when local orbitals are
Cr   present that have no corresponding smooth part.
Cr
Cr   The potential due to uniform background density rhobg is added
Cr   to v1 and v2 :   vbg = -(4pi/3)*(rmt^2-r^2)*rhobg
Cu Updates
Cu   10 Apr 19 Replace write statements with awrite ... buffers parallel
Cu   16 Jun 16 Addition of optional external potential vext
Cu   30 Sep 08 Return l=2 potential at nucleus in efg (A Svane)
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   19 Sep 02 added uniform bkg potential interactions (WRL)
Cu    8 Feb 02 rhoex and rhoec (T. Miyake)
Cu   14 Jan 02 rvexv and rvecv (T. Miyake)
Cu   15 Aug 01 Generates rvepsv and rvvxcv
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nr,nsp,nspc,lfoc,nlml,job,nbf
      double precision a,ceh,cofg,cofh,rfoc,rg,rhobg,rmt,z
      double precision bscal(2),bloc(nbf,3),vext(nr,nlml,nsp)
      real(8),target :: rho1(nr,nlml,nsp*nspc),rho2(nr,nlml,nsp*nspc)
      double precision rofi(nr),rwgt(nr),qmom(nlml),vval(nlml),gpotb(nlml),
     .  rhoc(nr,nsp*nspc),rhol1(nr,nlml,nsp),rhol2(nr,nlml,nsp),
     .  v1(nr,nlml,nsp*nspc),v1es(nr,nlml,nsp),
     .  v2(nr,nlml,nsp*nspc),v2es(nr,nlml,nsp)
      double precision efg(5)
      character*(*) buf
C ... For structures
!       include 'structures.h'
      type(str_atparms)::   s_atparms
C ... Dynamically allocated local arrays
      real(8),allocatable :: fl(:,:,:),yl(:,:),wp(:,:),wk(:,:,:),
     .  p(:,:),p2(:,:),r2(:),rhop(:,:,:),rhochs(:),rhonsm(:),rhocsm(:)
      complex(8),allocatable :: z1(:,:,:,:),z2(:,:,:,:)
      real(8),pointer :: rhoz1(:,:,:),rhoz2(:,:,:)
C ... Local parameters
      logical lvext
      integer stdo,ipr,ll,i,isp,ilm,l,lxcfun,nrml,isw,isw2,np,nspxc
      double precision afoc,ag,b,cof0,fac,gnu,pi,qcor1,qcor2,r,
     .  rhves1,rhves2,rvs1,rvs2,samh,sfac,smrhoc,srfpi,sum1,sum2,sumg,sumh,top,
     .  vcpn1,vcpn2,rhcvef1,ves1,vesc1,vesc2,vesn1,vesn2,vnucl,vsum,vtr,y0
      double precision a1(0:3),a2(0:3),qs(2),rep1c(2),rep1x(2),
     .  rep2c(2),rep2x(2),rvsm(2),rvtr(2),rvsx(2),rvtx(2),
     .  df(0:20),cof(nlml),xi(0:40),xx(2),rvext(2,3)
      procedure(integer) :: iprint,nglob
      procedure(real(8)) :: dot3,ddot,dlength

C --- Setup ---
C     print *, 'zero vval'
C     call dpzero(vval,nlml)
      call tcn('locpt2')
      stdo = nglob('stdo')
      ipr = iprint()
      isw = mod(job/1000,10)
      isw2 = mod(job/100000,10)*10
      lvext = mod(job/1000000,10) /= 0 ! Add external potential
      allocate(rhochs(nr*2),rhocsm(nr),rhonsm(nr))
      if (isw /= 0) then
        allocate (fl(nr,nlml,nsp+1))
      else
        allocate (fl(1,1,1))
      endif
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4*pi)
      y0 = 1d0/srfpi
      call stdfac(20,df)
      b = rmt/(dexp(a*nr-a)-1)
      ag = 1d0/rg
      afoc = 1d0/rfoc
      nrml = nr*nlml

C --- Make core and nucleus pseudodensities ---
      fac = 4*pi*(ag*ag/pi)**1.5d0
C ... Renormalize gaussian
      sumg = 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        sumg = sumg + rwgt(i)*gnu
      enddo
      if (dabs(sumg-1d0) > 1d-4) call awrit1('%N locpot (warning): '
     .  //'large gaussian, integral=%g',buf,80,stdo,sumg)
      sfac = 1d0/sumg
      fac = fac*sfac

C     sum1  = 0d0
C     qcor2 = 0d0
C     qcor1 = 0d0
      sumh = 0d0
C     Smooth nuc. and core rho, sm Hankel portion, true & smooth core q
      do  i = 1, nr
        r = rofi(i)
        gnu = fac * r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,afoc,xi,1)
        rhonsm(i) = -z*gnu
        smrhoc = srfpi*cofh*xi(0)*r*r
        rhocsm(i) = srfpi*cofg*gnu + smrhoc
        rhochs(i) = smrhoc
        sumh  = sumh + rwgt(i)*rhochs(i)
C       qcor1 = qcor1 + rwgt(i)*rhoc(i)
C       qcor2 = qcor2 + rwgt(i)*rhocsm(i)
C       sum1 = sum1 + rwgt(i)*rhonsm(i)
      enddo
      samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
      if (ipr >= 20 .and. dabs(samh) > 1d-6)
     .  call awrit3(' sm core charge = %;6d (sphere) + %;6d (spillout) = %;6d',
     .  buf,80,stdo,sumh,samh-sumh,samh)

C --- Noncollinear case: rotate rho1,rho2 into local quantization axis ---
      if (nspc == 2) then
C       call prrmsh('rho1 on entry',rofi,rho1,nr,nr,nlml*4)
C       call prrmsh('rho2 on entry',rofi,rho2,nr,nr,nlml*4)

        np = 122
        allocate(yl(np,nlml),wp(3,np))
        allocate(p(3,np),p2(np,3),r2(np))
        i = -np
        call fpiint(i,0,np,p,wp)  ! np
        call sanrg(.true.,-i,np,np,'locpot','np')
        call dmcpy(p,1,3,p2,np,1,np,3)
        call ropyln(np,p2,p2(1,2),p2(1,3),ll(nlml),np,yl,r2)
        deallocate (p,p2,r2)
        allocate(rhop(nr,np,nsp*nspc)) ! For density pointwise in sphere
        allocate(rhoz1(nr,nlml,nsp))   ! True valence density along local quantization axis
        allocate(rhoz2(nr,nlml,nsp))   ! Smooth density along local quantization axis
        allocate(z1(nr,np,nsp,nsp))    ! Eigenvectors of true density rotation
        allocate(z2(nr,np,nsp,nsp))    ! Eigenvectors of sm density rotation
        call fyl2p(nr,nlml,4,np,yl,rho1,rhop) ! Real space mesh
        call rotspa(4,nr*np,nr*np,nr*np,rhop,rhop,z1) ! Rotate to quantization axis
        z1(1,:,:,:) = z1(2,:,:,:)      ! rho(r=0) = 0: use rotation from 1st point
        call fp2yl(nr,nlml,nsp,np,wp,rhop,yl,0d0,rhoz1) ! Project to Ylm
        call fyl2p(nr,nlml,4,np,yl,rho2,rhop) ! Real space mesh
        call rotspa(4,nr*np,nr*np,nr*np,rhop,rhop,z2) ! Rotate to quantization axis
        call fp2yl(nr,nlml,nsp,np,wp,rhop,yl,0d0,rhoz2) ! Ylm repsn
        z2(1,:,:,:) = z1(2,:,:,:)      ! rho(r=0) = 0: use rotation from 1st point
        deallocate(rhop)
C       call dpzero(v1,nr*nlml*4)  ! not needed; rotspa(12,..) uses diagonal only
      else
        rhoz1 => rho1
        rhoz2 => rho2
      endif

C --- rhol1 = full electron density = rho1 + rhoc ---
      do  isp = 1, nsp
        call dpcopy(rhoz1(1,1,isp),rhol1(1,1,isp),1,nrml,1d0)
        call dpadd (rhol1(1,1,isp),rhoc(1,isp),1,nr,y0)
      enddo

C --- rhol2 = full smooth compensated density rho2+gval+gcor+gnuc ---
C     gval : qmom * compensating gaussians
C     rhonsm should integrate to -z exactly
C     Distribute core+nuclear charge equally over spins
      do  isp = 1, nsp
        do  ilm = 1, nlml
          l = ll(ilm)
          cof(ilm) = qmom(ilm)*4*pi/df(2*l+1)
          fac = sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
          do  i = 1, nr
            r = rofi(i)
            gnu = cof(ilm)*fac* r**(l+2) * dexp(-ag*ag*r*r)
            rhol2(i,ilm,isp) = rhoz2(i,ilm,isp) + gnu/nsp
          enddo
        enddo
        do  i = 1, nr
          rhol2(i,1,isp) = rhol2(i,1,isp) + y0/nsp*(rhocsm(i)+rhonsm(i))
        enddo
      enddo

C ... Combine separate spin densities for electrostatics
      call splrho(0,nsp,nr,nlml,rhoz1,rhoz2,rhoc)
      call splrho(20,nsp,nr,nlml,rhol1,rhol2,rhoc)

C ... Add background density to spherical rhol1 and rhol2
      do i = 1, nr
        rhol1(i,1,1) = rhol1(i,1,1) + srfpi*rhobg*rofi(i)**2
        rhol2(i,1,1) = rhol2(i,1,1) + srfpi*rhobg*rofi(i)**2
      enddo

C ... Sphere charges; also check sphere neutrality for safety
C     qv1+qcor1 = sum1+z  should be qv2+qcor2 = sum2+z
      s_atparms%qv1 = srfpi*ddot(nr,rwgt,1,rhoz1,1)
      s_atparms%qv2 = srfpi*ddot(nr,rwgt,1,rhoz2,1)
      a1 = 0 ; a2 = 0
      a1(0) = srfpi*ddot(nr,rwgt,1,rhoz1(1,1,nsp),1)
      a2(0) = srfpi*ddot(nr,rwgt,1,rhoz2(1,1,nsp),1)
      qcor1 = ddot(nr,rwgt,1,rhoc,1)
      qcor2 = ddot(nr,rwgt,1,rhocsm,1)
      s_atparms%qlocc = qcor1-qcor2
      if (nsp == 2) then
        s_atparms%alocc =  ddot(nr,rwgt,1,rhoc(1,2),1)
      endif
      s_atparms%qloc  = s_atparms%qv1 - s_atparms%qv2
      s_atparms%qcor1 = qcor1
      sum1  = srfpi*ddot(nr,rwgt,1,rhol1,1) - z
      sum2  = srfpi*ddot(nr,rwgt,1,rhol2,1) ! sm nuc charge is inside rhol2
      if (dabs(sum1-sum2) > 1d-6)
     .  call rx1('locpt2: sphere not neutral: charge = %d',sum1-sum2)
C     am = a1(0)-a2(0) ! integral mz (quantization axis z at each point)
C     bfam = double counting for external B field ...
      s_atparms%aloc = 0
      if (nsp == 1) then
        s_atparms%bfam = 0
C     Valid for nbf=1 and collinear case only!
      elseif (nspc == 1) then
        s_atparms%aloc(0) = a1(0)-a2(0)
        s_atparms%bfam = 0
        if (nbf > 0) s_atparms%bfam = a1(0)*bloc(1,3)
      else  !nspc == 2
        call magvec(nr,nlml,nsp,nspc,rwgt,rho1,a1)
        call magvec(nr,nlml,nsp,nspc,rwgt,rho2,a2)
        s_atparms%aloc(1:3) = a1(1:3) - a2(1:3)
        s_atparms%aloc(0) = dlength(3,s_atparms%aloc(1),1)
C       s_atparms%aloc(1:3) = s_atparms%aloc(1:3)/s_atparms%aloc(0); print *, s_atparms%aloc(1:3); stop
        if (nbf > 0) call rx('bfield not implemented in nc case')
      endif
      s_atparms%atrue = a1

C --- Solve Poisson equation for the true and smooth densities ---
C     call dpzero(vval, nlml)
C     Make v1 = Ves[rho1] = estat pot of true density rhol1 without nuclear contribution
C     Make v2 = Ves[rho2] = estat pot of sm density rhol2 which contains sm. nuc. rho
C     Note rvs1 = estat energy contains nuclear contribution int[rhol1 * (-2*Z/r)]
C     Boundary conditions are then:
C     v1(nr,ilm,1) = vval(ilm) + 2*z/rmt/y0 delta_ilm,1
C     v2(nr,ilm,1) = vval(ilm)
      allocate(wk(nr,nlml,nsp))
C     call poinsp(0d0,vval,nlml,a,b,v1,rofi,rhol1,wk,nr,rvs2,rhves1,vnucl,vsum)
      call poinsp(z,vval,nlml,a,b,v1,rofi,rhol1,wk,nr,rvs1,rhves1,vnucl,vsum)
C     Make v2 = Ves[rhol2 = rho2+gval+gcor+gnuc] (contains sm nuclear + core)
      call poinsp(0d0,vval,nlml,a,b,v2,rofi,rhol2,wk,nr,rvs2,rhves2,vnucl,vsum)
      deallocate(wk)

C     call prrmsh('rho1',rofi,rhol1,nr,nr,nlml)
C     call prrmsh('v1',rofi,v1,nr,nr,nlml)
C     call prrmsh('v2',rofi,v2,nr,nr,nlml)

C ... Electric field gradient at the nucleus
      efg = 0
      if (nlml >= 9 .and. z > 0.01d0) then
        do  ilm = 5, 9
          efg(ilm-4) = v1(5,ilm,1)/rofi(5)**2
        enddo
      endif

C --- gpotb = integrals of compensating gaussians times smooth ves + vext ---
C     vext is included in gpotb because rhol2 has n0~-n0 added to it (see above)
C     The corresponding energy term is sgpotb
C     sgpote computes separately multipole moments * vext
C     Term is not added to the energy because it is already contained in v2
C     In any event contributions made here should synchronize with those in smvxt
      s_atparms%sgpotb = 0d0; s_atparms%sgpote = 0d0
      do  ilm = 1, nlml
        l = ll(ilm)
        cof0 = 4*pi/df(2*l+1)
        fac = sfac*(ag*ag/pi)**1.5d0 * (2*ag*ag)**l
        sum1 = 0d0; sum2 = 0d0
        do  i = 1, nr
          r = rofi(i)
          gnu = cof0*fac* r**(l+2) * dexp(-ag*ag*r*r)
          sum1 = sum1 + rwgt(i)*v2(i,ilm,1)*gnu
          if (lvext) sum2 = sum2 + rwgt(i)*vext(i,ilm,1)*gnu
        enddo
        gpotb(ilm) = sum1 + sum2 ! Hartree + external, enters into ppi (see smvxt.f)
        s_atparms%sgpotb = s_atparms%sgpotb + qmom(ilm)*sum1
        s_atparms%sgpote = s_atparms%sgpote + qmom(ilm)*sum2
      enddo
C     call prmx('gpotb',gpotb,nlml,nlml,1)

C --- Contribution rho*vext to total energy, smooth and true parts ---
      rvext = 0
      if (lvext) then
C       Debugging check: copy v1 to vext and confirm energy = rvs1
C       call dcopy(nr*nlml*nsp,v1,1,vext,1)
C       do  np = 2, nr
C         vext(np,1,1:nsp) = vext(np,1,1:nsp) - 2d0*z/rofi(np)/y0
C       enddo
C       Debugging check: copy v2 to vext and confirm energy = rvs2
C       call dcopy(nr*nlml*nsp,v2,1,vext,1)
        do  isp = 1, nsp
          do  ilm = 1, nlml
            rvext(1,1) = rvext(1,1) + dot3(nr,rhol1(1,ilm,isp),vext(1,ilm,isp),rwgt(1))
            rvext(2,1) = rvext(2,1) + dot3(nr,rhol2(1,ilm,isp),vext(1,ilm,isp),rwgt(1))
          enddo
          rvext(1,2) = rvext(1,2) + y0*dot3(nr,rhoc,vext(1,1,isp),rwgt)
          rvext(2,2) = rvext(2,2) + y0*dot3(nr,rhocsm,vext(1,1,isp),rwgt)
          rvext(2,3) = rvext(2,3) - y0*dot3(nr,rhonsm,vext(1,1,isp),rwgt)
          rvext(1,3) = rvext(1,3) + z*vext(1,1,isp)*y0
        enddo
        rvext(2,1) = rvext(2,1) + rvext(2,3) - rvext(2,2)  ! subtract core, nuc
C       rvext(2,1) = rvext(2,1) - s_atparms%sgpote         ! subtract gaussians
        rvext(1,1) = rvext(1,1) - rvext(1,2)               ! subtract core

C        print *, 'q1,q2',ddot(nr,rwgt,1,rhol1,1)/y0,
C     .                   ddot(nr,rwgt,1,rhol2,1)/y0-ddot(nr,rwgt,1,rhonsm,1)

        if (ipr >= 50) then
          write(stdo,344)
  344     format(' Integrals of density with external potential:'/
     .      12x,'fvalence      core        n0~-n0         total        nucl')
          write(stdo,'(a,2f13.7,13x,2f13.7)') '   true',
     .      rvext(1,1), rvext(1,2), rvext(1,1)+rvext(1,2), rvext(1,3)
          write(stdo,'(a,5f13.7/)') ' smooth', rvext(2,1)-s_atparms%sgpote, rvext(2,2),
     .      s_atparms%sgpote, rvext(2,1)+rvext(2,2), rvext(2,3)
        endif
        s_atparms%rvext(1:2,1:2) = rvext(1:2,1:2)
      endif

C --- Electrostatic integrals involving spherical terms only ---
C     NB: v1 and v2 contain external potential, if it exists
C         v2 includes nuclear part of rho, v1 does not
      vesc1 = 0d0  ! int (true core charge) * v1
      vesn2 = 0d0  ! int (sm   nuc charge)  * v2
      vesc2 = 0d0  ! int (sm   core charge) * v2
      vnucl = 0d0  ! int (true  val charge) * (1/R-1/r) (BC=0 for vnucl)
      do  i = 2, nr
        ves1 = y0*v1(i,1,1) - 2*z/rofi(i)
        vesc1 = vesc1 + rwgt(i)*rhoc(i,1)*ves1
        vesn2 = vesn2 + rwgt(i)*rhonsm(i)*y0*v2(i,1,1)
        vesc2 = vesc2 + rwgt(i)*rhocsm(i)*y0*v2(i,1,1)
        vnucl = vnucl + rwgt(i)*rhol1(i,1,1)*(2d0/rofi(i)-2d0/rmt)/y0
      enddo
C     Debugging check make v(r) ... show that vnucl = v(1) if v(rmt) = 0
C     call lscoul(rofi,nr,nsp,1,0d0,0d0,rhol1,v1)
C     cof0 = 2/y0*ddot(nr,rwgt,1,rhol1,1)/rofi(nr) ! 2Q/rmt
C     print *, 'v(rmt)',v1(nr,1,1)*y0 - cof0
C     print *, 'v(1)',  v1(1,1,1)*y0 - cof0, vnucl
C     stop

C     vnucl from preceding = ves(0) from rhol1, corresponding to
C     B.C. ves(rmt)=0.  Proper BC is (cf call to poinsp above):
C     ves(rmt) = v1(nr,1,1)*y0 = vval(1)*y0 + 2*z/rmt
      vnucl = vnucl + 2*z/rmt + y0*vval(1)
      vesn1 = -z * vnucl

C ... Valence density times electrostatic potential
      s_atparms%vales1 = rvs1-vesc1            ! Remove rhoc * ves to get valence only
      s_atparms%vales2 = rvs2-vesn2-vesc2      ! Remove core here also, and nuclear part
      s_atparms%vvesat = s_atparms%vales1-s_atparms%vales2

C ... Core plus nucleus times estatic potential
      vcpn1  = vesc1 + vesn1
      vcpn2  = vesn2 + vesc2
      s_atparms%cpnves = vcpn1 - vcpn2

C ... Subtract background before doing exchange correlation
      do i = 1, nr
        rhol1(i,1,1) = rhol1(i,1,1)-srfpi*rhobg*rofi(i)**2
        rhol2(i,1,1) = rhol2(i,1,1)-srfpi*rhobg*rofi(i)**2
      enddo

C ... Restore separate spin densities; copy estat pot to spin2
      if (nsp == 2) then
        call splrho(1,nsp,nr,nlml,rhoz1,rhoz2,rhoc)
        call splrho(21,nsp,nr,nlml,rhol1,rhol2,rhoc)
        call dcopy(nrml,v1,1,v1(1,1,2),1)
        call dcopy(nrml,v2,1,v2(1,1,2),1)
      endif

C ... Keep Hartree part
      call dcopy(nrml*nsp,v1,1,v1es,1)
      call dcopy(nrml*nsp,v2,1,v2es,1)

C     Case make vxc for spin averaged density:
C     Combine densities, treat as though nsp=1
      nspxc = nsp
      if (nsp == 2 .and. mod(job/100,10) == 2) then
        call splrho(0,nsp,nr,nlml,rhoz1,rhoz2,rhoc)
        call splrho(20,nsp,nr,nlml,rhol1,rhol2,rhoc)
        nspxc = 1
      endif

C ... Generate valence-only rvvxcv,rvepsv,rvexv,rvecv  (uses v1 as work array)
C     call pshpr(ipr-31)
      if (isw == 1) then
        if (nspxc /= nsp) call rx('locpt2 not ready for job=200')
        call pshpr(max(ipr-31,min(ipr,10)))
        lxcfun = nglob('lxcf')
        call vxcnsp(isw2,rofi,nr,rwgt,nlml,nspxc,rhoz1,0,lxcfun,bscal(1),0,
     .    0,0,0,0,s_atparms%rep1,rep1x,rep1c,s_atparms%rmu1,v1,fl,qs,buf)
C       call prrmsh('rho2 before',rofi,rhoz2,nr,nr,nlml)
        call vxcnsp(isw2,rofi,nr,rwgt,nlml,nspxc,rhoz2,0,lxcfun,bscal(2),0,
     .    0,0,0,0,s_atparms%rep2,rep2x,rep2c,s_atparms%rmu2,v1,fl,qs,buf)
C       call prrmsh('rho2 after',rofi,rhoz2,nr,nr,nlml)
        s_atparms%rvvxcv = s_atparms%rmu1(1) - s_atparms%rmu2(1)
     .                   + s_atparms%rmu1(2) - s_atparms%rmu2(2)
        s_atparms%rvepsv = s_atparms%rep1(1) - s_atparms%rep2(1)
     .                   + s_atparms%rep1(2) - s_atparms%rep2(2)
        s_atparms%rvexv  = rep1x(1) - rep2x(1) + rep1x(2) - rep2x(2)
        s_atparms%rvecv  = rep1c(1) - rep2c(1) + rep1c(2) - rep2c(2)
        call dcopy(nrml*nspxc,v1es,1,v1,1)
        call poppr
      endif

C --- Add xc potentials to v1 and v2 ---
C     call pshpr(ipr-11)
      call pshpr(max(ipr-11,min(ipr,10)))
      lxcfun = nglob('lxcf')
      call dpzero(s_atparms%focexc,2)
      call dpzero(s_atparms%focex,2)
      call dpzero(s_atparms%focec,2)
      call dpzero(s_atparms%focvxc,2)
      if (iprint() >= 30) call awrit0(' Exchange for true density:',buf,100,stdo)
C     call prrmsh('rho1 before',rofi,rhol1,nr,nr,nlml)
      call vxcnsp(isw+isw2,rofi,nr,rwgt,nlml,nspxc,rhol1,0,lxcfun,
     .  bscal(1),0,0,0,0,0,s_atparms%rep1,rep1x,rep1c,s_atparms%rmu1,v1,fl,qs,buf)
C     call prrmsh('rho1 after',rofi,rhol1,nr,nr,nlml)
      if (isw == 1) then
       call dpzero(xx,2)
       call vxcns5(1,31,'rhov*vxc',nlml,nspxc,nr,rofi,rwgt,rhoz1,fl,xi,xx,buf)
       s_atparms%rvvxc = xx(1) + xx(2)
       do  isp = 1, nspxc
         call vxcns5(1,31,'rhov*exc',nlml,1,nr,rofi,rwgt,rhoz1(1,1,isp),
     .     fl(1,1,3),xi,xx(isp),buf)
       enddo
       s_atparms%rveps = xx(1) + xx(2)
      endif

C     call prrmsh('v1',rofi,v1,nr,nr,nlml)
C     If no core treatment v2 += vxc(rho2)
      if (lfoc == 0) then
        call vxcnsp(isw+isw2,rofi,nr,rwgt,nlml,nspxc,rhoz2,0,lxcfun,
     .    bscal(2),0,0,0,0,0,s_atparms%rep2,rep2x,rep2c,s_atparms%rmu2,v2,fl,qs,buf)
C     Otherwise v2 += vxc(rho2 + sm core), directly or perturbatively:
      else if (lfoc == 1) then
        if (iprint() >= 30) call awrit0(' Exchange for smooth density, foca=1:',buf,100,stdo)
        do  isp = 1, nspxc
          call dpadd(rhoz2(1,1,isp),rhochs,1,nr,y0/nspxc)
        enddo
        call vxcnsp(isw+isw2,rofi,nr,rwgt,nlml,nspxc,rhoz2,0,lxcfun,
     .    bscal(2),0,0,0,0,0,s_atparms%rep2,rep2x,rep2c,s_atparms%rmu2,v2,fl,qs,buf)
        do  isp = 1, nspxc
          call dpadd(rhoz2(1,1,isp),rhochs,1,nr,-y0/nspxc)
        enddo
C     v2 += vxc(rho2) + dv; focvxc= int (rho2*dv); dv=dvxc/drho * rhochs
      else if (lfoc == 2) then
        if (iprint() >= 30) call awrit0(' Exchange for smooth density, foca=2:',buf,100,stdo)
        if (nspxc == 2) then
          call dcopy(nr,rhochs,1,rhochs(1+nr),1)
        endif
C       call prrmsh('rho2',rofi,rhoz2,nr,nr,nlml*nspxc)
        call vxcnsp(isw+isw2,rofi,nr,rwgt,nlml,nspxc,rhoz2,1,lxcfun,bscal(2),rhochs,s_atparms%focexc,
     .    s_atparms%focex,s_atparms%focec,s_atparms%focvxc,s_atparms%rep2,rep2x,rep2c,s_atparms%rmu2,v2,fl,qs,buf)
        if (ipr >= 50) write(stdo,941) sumh,
     .    s_atparms%focexc(1)+s_atparms%focexc(2),s_atparms%focvxc(1)+s_atparms%focvxc(2)
  941   format(' smH core xc integrals',3f12.6)
      else
        call rxi('locpt2: cannot handle lfoc = ',lfoc)
      endif

      if (isw == 1) then
       call vxcns5(1,31,'rhov*vxc',nlml,nspxc,nr,rofi,rwgt,rhoz2,fl,xi,xx,buf)
       s_atparms%rvvxc = s_atparms%rvvxc - (xx(1) + xx(2))
       do  isp = 1, nspxc
         call vxcns5(1,31,'rhov*exc',nlml,1,nr,rofi,rwgt,rhoz2(1,1,isp),
     .     fl(1,1,3),xi,xx(isp),buf)
       enddo
       s_atparms%rveps = s_atparms%rveps - (xx(1) + xx(2))
      endif
      call poppr
C     call prrmsh('v2',rofi,v2,nr,nr,nlml)

      if (nsp == 2 .and. mod(job/100,10) == 2) then
        call splrho(1,nsp,nr,nlml,rhoz1,rhoz2,rhoc)
        call splrho(21,nsp,nr,nlml,rhol1,rhol2,rhoc)
        call dcopy(nrml,v1,1,v1(1,1,2),1)
        call dcopy(nrml,v2,1,v2(1,1,2),1)
        s_atparms%rep1(2) = s_atparms%rep1(1)/2
        s_atparms%rep1(1) = s_atparms%rep1(2)
        rep1x(2) = rep1x(1)/2; rep1x(1) = rep1x(2)
        rep1c(2) = rep1c(1)/2; rep1c(1) = rep1c(2)
        s_atparms%rmu1(2) = s_atparms%rmu1(1)/2
        s_atparms%rmu1(1) = s_atparms%rmu1(2)
        s_atparms%rep2(2) = s_atparms%rep2(1)/2
        s_atparms%rep2(1) = s_atparms%rep2(2)
        rep2x(2) = rep2x(1)/2; rep2x(1) = rep2x(2)
        rep2c(2) = rep2c(1)/2; rep2c(1) = rep2c(2)
        s_atparms%rmu2(2) = s_atparms%rmu2(1)/2
        s_atparms%rmu2(1) = s_atparms%rmu2(2)
      endif

C --- Integrals over core times effective potential ---
      rhcvef1 = 0d0
      do  isp = 1, nsp
      do  i = 2, nr
        ves1 = y0*v1(i,1,isp) - 2*z/rofi(i)
        rhcvef1 = rhcvef1 + rwgt(i)*rhoc(i,isp)*ves1
      enddo
      enddo

C --- Integrals involving the full nonspherical estat + xc potential ---
      s_atparms%rhovv1 = 0d0; s_atparms%rhovv2 = 0d0
      s_atparms%rhvxt1 = 0d0; s_atparms%rhvxt2 = 0d0
C      if (ipr >= 45 .and. nsp == 1) write(stdo,351)
C      if (ipr >= 45 .and. nsp == 2) write(stdo,353)

      if (ipr >= 45 .and. nsp == 1)
     .  call awrit0('%N ilm%9frho*vtrue%7frho*vsm',buf,100,stdo)
      if (ipr >= 45 .and. nsp == 2) call awrit0(
     .  '%N ilm%19frho*vtrue%30frho*vsm%N%13fspin1%7fspin2%7ftot%11fspin1%7fspin2%7ftot',
     .  buf,160,stdo)

      do  ilm = 1, nlml
        rvtr = 0d0; rvsm = 0d0; rvtx = 0d0; rvsx = 0d0
        do  isp = 1, nsp
          do  i = 2, nr
            vtr = v1(i,ilm,isp)
            if (ilm == 1) vtr = vtr - srfpi*2*z/rofi(i)
            rvtr(isp) = rvtr(isp) + rwgt(i)*rhoz1(i,ilm,isp)*vtr
            rvsm(isp) = rvsm(isp) + rwgt(i)*rhoz2(i,ilm,isp)*v2(i,ilm,isp)
            rvtx(isp) = rvtx(isp) + rwgt(i)*rhoz1(i,ilm,isp)*vext(i,ilm,isp)
            rvsx(isp) = rvsx(isp) + rwgt(i)*rhoz2(i,ilm,isp)*vext(i,ilm,isp)
          enddo
          s_atparms%rhovv1 = s_atparms%rhovv1 + rvtr(isp)
          s_atparms%rhovv2 = s_atparms%rhovv2 + rvsm(isp)
          s_atparms%rhvxt1 = s_atparms%rhvxt1 + rvtx(isp)
          s_atparms%rhvxt2 = s_atparms%rhvxt2 + rvsx(isp)
        enddo
        top = dmax1(dabs(rvsm(1)),dabs(rvtr(1)))
        if (ipr < 45) cycle
        if (top >= 1d-6 .and. nsp == 1) then
C          write(stdo,350) ilm,rvtr(1),rvsm(1)
C          if (lvext) write(stdo,354) rvtx(1),rvsx(1)
          call awrit3('%,4i   %;15,6D%;15,6D',buf,100,stdo,ilm,rvtr(1),rvsm(1))
          if (lvext) call awrit2(' r*vext%;15,6D%;15,6D',
     .      buf,100,stdo,rvtx(1),rvsx(1))
        endif
        if (top >= 1d-6 .and. nsp == 2) then
C          write(stdo,352) ilm,rvtr(1),rvtr(2),rvtr(1)+rvtr(2),
C     .      rvsm(1),rvsm(2),rvsm(1)+rvsm(2)
          call awrit7('%,4i   %;12,6D%;12,6D%;12,6D  %;12,6D%;12,6D%;12,6D',
     .      buf,100,stdo,ilm,rvtr(1),rvtr(2),rvtr(1)+rvtr(2),rvsm(1),rvsm(2),rvsm(1)+rvsm(2))
C          if (lvext) write(stdo,355) rvtx(1),rvtx(2),rvtx(1)+rvtx(2),
C     .      rvsx(1),rvsx(2),rvsx(1)+rvsx(2)
          if (lvext) call awrit7(' r*vext%;12,6D%;12,6D%;12,6D  %;12,6D%;12,6D%;12,6D',
     .      buf,100,stdo,ilm,rvtx(1),rvtx(2),rvtx(1)+rvtx(2),rvsx(1),rvsx(2),rvsx(1)+rvsx(2))

        endif
C  350   format(i4,3x,2f15.6)
C  352   format(i4,3x,3f12.6,2x,3f12.6,2x)
C  351   format(/' ilm',09x,'rho*vtrue',07x,'rho*vsm')
C  353   format(/' ilm',19x,'rho*vtrue',30x,'rho*vsm'/13x,
C     .    'spin1',7x,'spin2',7x,'tot',11x,'spin1',7x,'spin2',7x,'tot')
C  354   format(' r*vext',2f15.6)
C  355   format(' r*vext',3f12.6,2x,3f12.6,2x)
      enddo

C ... Smooth xc potential includes foca head; undo in integral
      s_atparms%rhovv2  = s_atparms%rhovv2 - s_atparms%focvxc(1)-s_atparms%focvxc(2)
      s_atparms%rhoexc(1) = s_atparms%rep1(1) - s_atparms%rep2(1)
      s_atparms%rhoexc(2) = s_atparms%rep1(2) - s_atparms%rep2(2)
      s_atparms%rhoex(1)  = rep1x(1) - rep2x(1)
      s_atparms%rhoex(2)  = rep1x(2) - rep2x(2)
      s_atparms%rhoec(1)  = rep1c(1) - rep2c(1)
      s_atparms%rhoec(2)  = rep1c(2) - rep2c(2)
      s_atparms%rhovxc(1) = s_atparms%rmu1(1) - s_atparms%rmu2(1)
      s_atparms%rhovxc(2) = s_atparms%rmu1(2) - s_atparms%rmu2(2)
      s_atparms%rhcvef1  = rhcvef1
      s_atparms%rhovv2 = s_atparms%rhovv2 + s_atparms%sgpotb
      s_atparms%rhvxt2 = s_atparms%rhvxt2 + s_atparms%sgpote  ! Should no longer add

C      call info5(50,0,0,' locpot rhovv1 = %,6;6d rhovv2 = %,6;6d'//
C     .  ' s_atparms%sgpotb = %,6;6d rhovef = %,6;6d',s_atparms%rhovv1,s_atparms%rhovv2,s_atparms%sgpotb,s_atparms%rhovv1-s_atparms%rhovv2,0)

C --- Charges, printout ---
      if (ipr >= 45) then
C        write(stdo,251)
C        write(stdo,250) s_atparms%rep1(1)+s_atparms%rep1(2),s_atparms%rep2(1)+s_atparms%rep2(2),
C     .                  s_atparms%rhoexc(1)+s_atparms%rhoexc(2),
C     .                  s_atparms%rmu1(1),s_atparms%rmu2(1),s_atparms%rhovxc(1)
C        if (nsp == 2) write(stdo,253) s_atparms%rmu1(2),s_atparms%rmu2(2),s_atparms%rhovxc(2),
C     .                  s_atparms%rmu1(1)+s_atparms%rmu1(2),s_atparms%rmu2(1)+s_atparms%rmu2(2),
C     .                  s_atparms%rhovxc(1)+s_atparms%rhovxc(2)
C        write(stdo,252) s_atparms%vales1,s_atparms%vales2,s_atparms%vales1-s_atparms%vales2,
C     .    s_atparms%rhovv1,s_atparms%rhovv2,s_atparms%rhovv1-s_atparms%rhovv2,s_atparms%qv1,
C     .    s_atparms%qv2,s_atparms%qloc
C        if (nsp == 2) write(stdo,254) a1(0),a2(0),a1(0)-a2(0),s_atparms%alocc
C        if (nspc == 2) then
C          write(stdo,255) '   Mtrue: ', s_atparms%atrue(1:3),  '  |Mtrue|',s_atparms%atrue(0)
C          write(stdo,256) '   Mloc:  ', s_atparms%aloc(1:3) ! ,'  M-int mz',s_atparms%aloc(0)-am
C        endif
C        write(stdo,257) qcor1,qcor2,s_atparms%qlocc
CC       call info2(50,0,0,' (n0~-n0) Ves%;12,5D',s_atparms%sgpotb,0)
C
C  251   format(/' local terms:     true',11x,'smooth',9x,'local')
C  250   format(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6)
C  253   format(' spin2:   ',3f15.6/' total:   ',3f15.6)
C  252   format(' val*ves  ',3f15.6/' val*vef  ',3f15.6/' val chg: ',3f15.6)
C  254   format(' val mom: ',3f15.6,'    core:',f11.6)
C  255   format(a,           3f15.6,a,f11.6)
C  256   format(a,           3f15.6,a,f10.6)
C  257   format(' core chg:',3f15.6)

        call awrit0('%N local terms:     true%11fsmooth%9flocal',buf,160,stdo)
        call awrit6(' rhoeps:  %;15,6D%;15,6D%;15,6D%N rhomu:   %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%rep1(1)+s_atparms%rep1(2),s_atparms%rep2(1)+s_atparms%rep2(2),
     .    s_atparms%rhoexc(1)+s_atparms%rhoexc(2),
     .    s_atparms%rmu1(1),s_atparms%rmu2(1),s_atparms%rhovxc(1))
        if (nsp == 2)
     .    call awrit6(' spin2:   %;15,6D%;15,6D%;15,6D%N total:   %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%rmu1(2),s_atparms%rmu2(2),s_atparms%rhovxc(2),
     .    s_atparms%rmu1(1)+s_atparms%rmu1(2),s_atparms%rmu2(1)+s_atparms%rmu2(2),
     .    s_atparms%rhovxc(1)+s_atparms%rhovxc(2))
        call awrit6(' val*ves  %;15,6D%;15,6D%;15,6D%N val*vef  %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%vales1,s_atparms%vales2,s_atparms%vales1-s_atparms%vales2,
     .    s_atparms%rhovv1,s_atparms%rhovv2,s_atparms%rhovv1-s_atparms%rhovv2)
        call awrit3(' val chg: %;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    s_atparms%qv1,s_atparms%qv2,s_atparms%qloc)
        if (nsp == 2) call awrit4(' val mom: %;15,6D%;15,6D%;15,6D    core:%;11,6D',buf,160,stdo,
     .    a1(0),a2(0),a1(0)-a2(0),s_atparms%alocc)
        if (nspc == 2) then
          call awrit2('   Mtrue: %3;15,6D  |Mtrue|%;11,6D',buf,160,stdo,
     .      s_atparms%atrue(1:3),s_atparms%atrue(0))
          call awrit1('   Mloc:  %3;15,6D',buf,160,stdo,s_atparms%aloc(1:3))
        endif
        call awrit3(' core chg:%;15,6D%;15,6D%;15,6D',buf,160,stdo,
     .    qcor1,qcor2,s_atparms%qlocc)
      endif

      deallocate (fl,rhochs,rhocsm,rhonsm)

      if (nspc == 2) then

        allocate(rhop(nr,np,nsp*nspc)) ! For density pointwise in sphere

C       Restore rho1 if vxcnsp may have altered it
        if (isw2 == 30) then
        call fyl2p(nr,nlml,4,np,yl,rhoz1,rhop) ! rhoz1 -> real space mesh
        call rotspa(12,nr*np,nr*np,nr*np,rhop,rhop,z1) ! Rotate to global axis
        call fp2yl(nr,nlml,4,np,wp,rhop,yl,0d0,rho1) ! Project to Ylm
C       call prrmsh('rho1 on exit',rofi,rho1,nr,nr,nlml*4)
        endif

C       Rotate v1 to global axis
        call fyl2p(nr,nlml,4,np,yl,v1,rhop) ! v1 -> real space mesh
        call rotspa(12,nr*np,nr*np,nr*np,rhop,rhop,z1) ! Rotate to global axis
        call fp2yl(nr,nlml,4,np,wp,rhop,yl,0d0,v1) ! Project to Ylm
C       call prrmsh('v1 on exit',rofi,v1,nr,nr,nlml*4)

C       Restore rho1 to to global axis if vxcnsp may have altered it
        if (isw2 == 30) then
        call fyl2p(nr,nlml,4,np,yl,rhoz2,rhop) ! rhoz2 -> real space mesh
        call rotspa(12,nr*np,nr*np,nr*np,rhop,rhop,z2) ! Rotate to global axis
        call fp2yl(nr,nlml,4,np,wp,rhop,yl,0d0,rho2) ! Project to Ylm
C       call prrmsh('rho2 on exit',rofi,rho2,nr,nr,nlml*4)
        endif

C       Rotate v2 to global axis
        call fyl2p(nr,nlml,4,np,yl,v2,rhop) ! v2 -> real space mesh
        call rotspa(12,nr*np,nr*np,nr*np,rhop,rhop,z2) ! Rotate to global axis
        call fp2yl(nr,nlml,4,np,wp,rhop,yl,0d0,v2) ! Project to Ylm
C       call prrmsh('v2 on exit',rofi,v2,nr,nr,nlml*4)

        deallocate(yl,wp,rhop,rhoz1,rhoz2,z1,z2)

      endif

C ... Add external potential to v1 and v2
      if (lvext) then
        call daxpy(nr*nlml*nsp,1d0,vext,1,v1,1)
        call daxpy(nr*nlml*nsp,1d0,vext,1,v2,1)
      endif

      call tcx('locpt2')
      end
