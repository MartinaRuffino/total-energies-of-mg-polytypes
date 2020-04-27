C#define F90
      subroutine vxcnsp(isw,ri,nr,rwgt,nlm,nsp,rl,lpert,lxcfun,bscal,rc,
     .  focexc,focex,focec,focvxc,reps,repsx,repsc,rmu,vl,fl,qs,buf)
C- Add vxc to potential in sphere and make integrals rho*vxc,rho*exc
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1s digit
Ci         : 1, make fl, otherwise not
Ci         :10s digit
Ci         : 0 calculates LDA XC energy, potential for density as is.
Ci         : 1 reset points of negative density to positive density
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci         : 3 reset points of negative density to positive density,
Ci         :   and alter rl if any rho<0 encountered
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rl    :full charge density * r**2
Ci   lpert :0 normal treatment of core XC potential, included with valence
Ci         :1 perturbative treatment of core XC potential
Ci   lxcfun:defines xc functional; see lxcf in structures.h
Ci   rc    :l=0 charge density*4*pi*ri**2 to be added in pert. theory
Ci         :Not used unless 10000s digit of lxcfun is set.
Ci   buf   :Buffer to contain redirected standard output.
Ci         :This routine uses awrite to write output to stdout.
Ci         :It gets redirected to buf if buffer mode is turned on
Co Outputs
Co   focex :(100s digit lxcfun): integral rc * vx
Co   focec :(100s digit lxcfun): integral rc * vc
Co   focexc:(100s digit lxcfun): integral rc * vxc
Co   focvxc:(100s digit lxcfun): integral rc * (dvxc/drho * rho)
Co         : for these integrals, see Remarks
Co   reps  :integral rl*exc
Co   repsx :integral rl*ex
Co   repsc :integral rl*ec
Co   rmu   :integral rl*vxc
Co   vl    :local exchange-correlation potential added to vl
Co   fl    :fl(:,:,1:nsp) = vxc by L
Co         :fl(:,:,1+nsp) = exc by L
Co         :fl is made only if isw is nonzero
Co   qs    :spin-resolved charge
Cl Local variables
Cl   yl    :YL(rp)
Cl   lxcfl :If nonzero, evaluate a local XC functional through vxcnsl
Cl         :If >1000, the libxc library is used.
Cl   lxcfnl:If nonzero, indicates a GGA from from the libxc library
Cl         :is to be evaluated through vxcnls
Cl   lxcg  :If nonzero, flags that this functional requires gradients
Cl         :If not a libxc functional, lxcg specifies which GGA
Cr Remarks
Cr   Perturbation treatment:
Cr     vxc[rho + rc] = vxc[rho] + rc * (dvxc/drho)
Cr   Correction to int rho * vxc:
Cr     int rho * vxc = int rho*vxc[rho] + int rho*rc*(dvxc/drho)
Cr     Correction = focvxc = int rho * rc * (dvxc/drho)
Cr   Correction to int rho * exc:
Cr     Corr = int rho * rc * (dexc/drho) = int rho * rc * (vxc-exc)/rho
Cr          = int rc * (vxc-exc) = focexc - int rc exc
Cr     Second term is not retained.
Cu Updates
Cu   10 Apr 19 Replace info or write statements with awrite ... buffers parallel
Cu   02 Jul 15 Remove fixed dimensioning nlmx
Cu   31 Jul 14 XC B-field can be scaled by bscal
Cu   09 Dec 13 First cut at calling libxc functional for XC potential
Cu   01 Dec 13 Enable calls to libxc library.  LDA only for now.
Cu   06 Apr 09 Routine calls vxcnls to handle some GGA potentials
Cu   14 Mar 04 Optionally Makes fl.  New argument list
Cu   14 Jun 02 repsx and repsc (T. Miyake)
Cu    8 Feb 02 focex and focec (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    3 Jun 00 Adapted from old FP vxcnsp; pert. corr from nfp vxcdnsp.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer isw,nr,nlm,lpert,lxcfun,nsp
      double precision bscal,ri(nr),reps(2),rmu(2),repsx(2),repsc(2),
     .  rwgt(nr),rc(nr),rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,1+nsp)
      double precision focexc(2),focex(2),focec(2),focvxc(2)
      double precision qs(2)
      character*(*) buf
C ... Dynamically allocated local arrays
      real(8), allocatable :: yl(:,:),ylwp(:,:),gyl(:,:,:)
C ... Local parameters
      integer nnn
      parameter(nnn=426)
      integer stdo,ipr,lmax,np,nph,nth,nxl(0:7)
      integer xctype(2),lxcfl,lxcg,lxcfnl
      double precision p(3,nnn),wp(nnn),p2(nnn*3),r2(nnn),fpi
      procedure(integer) :: lgunit,ll

      data nxl /-12,-20,-32,-32,-62,-92,-122,-122/

      stdo = lgunit(1)
      call getpr(ipr)
      fpi = 16d0*datan(1d0)
      if (lxcfun >= 1000) then   ! a libxc functional.  Check if need gradients
        call xc_libxc_type(lxcfun,0,xctype)
        if (xctype(1) > 2 .or. xctype(2) > 2)
     .    call rx('vxcnsp: functional not implemented')
        if (xctype(1) == 2 .or. xctype(2) == 2) then ! This is GGA
          lxcfl = 0        ! No potential calculated by vxcnsl (LDA)
          lxcfnl = lxcfun  ! Potential calculated by vxcnls (GGA)
          lxcg = 1         ! Require gradients
        else ! This is LDA
          lxcfl = lxcfun ! Evaluate XC functional through vxcnsl (LDA)
          lxcfnl = 0     ! No potential calculated by vxcnsl (GGA)
          lxcg = 0       ! No gradients
        endif
      else
        lxcfl = mod(lxcfun,100) ! Evaluate XC functional through vxcnsl
        lxcfnl = 0              ! No libxc potential calculated by vxcnls
        lxcg = mod(lxcfun/100,100) ! nonlocal part of GGA by vxcnls
      endif

C ... Create angular integration mesh
      lmax = ll(nlm)
C     if (lxcg /= 0) lmax = lmax+1
      if (lmax > 6) then
        nth = 2*lmax+2
        nph = 0
      else
        nth = nxl(lmax)
        nph = 0
      endif
      call fpiint(nth,nph,np,p,wp)
C      if (ipr >= 30) write (stdo,1) nth,nph,np,nr
C    1 format(' mesh:   nth,nph=',2I4,'   gives',i4,'  angular points,','   nrad=',i4)
      if (ipr >= 30) call awrit4(' mesh:   nth,nph= %i %i   gives %i angular points, %i radial',
     .  buf,100,stdo,nth,nph,np,nr)
      if (np > nnn) call rxi('vxcnsp: increase nnn, need',np)
C     if ((lmax+2)**2*np > nlmx*nnn) call rx('vxcnsp: increase nlm')

      rmu(2)  = 0
      reps(2) = 0
      repsx(2) = 0
      repsc(2) = 0

C ... Scale rl to true density
      call vxcns3(nr,nlm,nsp,ri,rl,1)
      if (lpert /= 0) then
        call vxcns3(nr,1,1,ri,rc,1)
        call dscal(nr,1/fpi,rc,1)
      endif

C ... Spherical harmonics (unit radius) ... to be overwritten by Yl * wp
      call dmcpy(p,1,3,p2,np,1,np,3)
C     GGA case:  Need both yl*wp and YL to lmax+1 to make grad YL
      if (lxcg /= 0) then
        allocate (yl(np,(lmax+2)**2),ylwp(np,(lmax+2)**2))
        call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax+1,np,yl,r2)
        call dcopy(size(yl),yl,1,ylwp,1)
      else
        allocate (ylwp(np,nlm))
        call ropyln(np,p2,p2(1+np),p2(1+2*np),lmax,np,ylwp,r2)
      endif

C ... Set up the density (checking for negative values)
C     Generate the local exchange-correlation potential
C     Return ylwp = YL * wp
      call vxcnsl(isw,ri,nr,rwgt,np,wp,ylwp,nlm,nsp,rl,rc,lxcfl,bscal,
     .  lpert,focexc,focex,focec,focvxc,reps,repsx,repsc,rmu,
     .  vl,fl,qs,buf)

C --- GGA: correction to LD potential, or substitute for LDA ---
      if (lxcg /= 0) then

C       Gradient of spherical harmonics
        allocate (gyl(np,nlm,3))
        call ropylg(1,lmax,nlm,np,np,p2,p2(1+np),p2(1+2*np),r2,yl,gyl)
C        allocate(grp(nr*np*nsp,3),agrp(nr*np*nsp),ggrp(nr*np*nsp))
C        allocate(rp(nr,np,nsp))
C        allocate(agrl(nr*nlm*nsp),vxcnl(nr*nlm*nsp),excnl(nr*nlm*nsp))
        call vxcnls(ri,0,nr,np,nlm,nsp,yl,gyl,ylwp,rwgt,wp,rl,
     .    lxcfun,bscal,vl,reps,rmu,buf)
        deallocate(yl,gyl)

      endif
      deallocate(ylwp)

C     Debugging
C      print *, bscal
C      print *, vl(1,1,1)+vl(1,1,2),vl(1,1,1)-vl(1,1,2)
C      print *, vl(nr,1,1)+vl(nr,1,2),vl(nr,1,1)-vl(nr,1,2)
C      stop

C ... Undo the rl scaling
      call vxcns3(nr,nlm,nsp,ri,rl,0)
      if (lpert /= 0) then
        call vxcns3(nr,1,1,ri,rc,0)
        call dscal(nr,fpi,rc,1)
      endif

      end

      subroutine vxcnsl(isw,ri,nr,rwgt,np,wp,yl,nlm,nsp,rl,rc,lxcf,bscal,
     .  lpert,focexc,focex,focec,focvxc,rep,repx,repc,rmu,vl,fl,qs,buf)
C- Make vxc, rep and rmu in sphere, for local XC functionals.
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :1s digit
Ci         : 1, make fl, otherwise not
Ci         :10s digit
Ci         : 0 calculated LDA for density as is.
Ci         : 1 reset points of negative density to positive density
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci         : 3 reset points of negative density to positive density,
Ci         :   and alter rl if any rho<0 encountered
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   np    :number of spherical mesh points
Ci   wp    :spherical mesh weights
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   rl    :charge density.
Ci         :Note: rl will be altered if points of negative density
Ci         :and 10s digit of isw is 3.
Ci   rc    :l=0 charge density*4*pi*ri**2 to be added in pert. theory
Ci   lxcf  :specifies local exchange-correlation functional
Ci         :  0    No LDA, no potential generated.
Ci         :       Reset negative spots of density and return yl
Ci         :  1    Ceperly-Alder
Ci         :  2    Barth-Hedin (ASW fit)
Ci         :  3,4  LD part of PW91 and PBE
Ci   lpert :Make perturbation integrals focexc and focvxc
Cio Inputs/Outputs
Cio  yl    :Spherical harmonics YL(rp)
Cio        :yl is OVERWRITTEN by YL * wp on output
Co Outputs
Co   qs    :spin-resolved charge
Co   focexc:(lpert only): integral rc * vxc
Co   focex :(lpert only): integral rc * vx
Co   focec :(lpert only): integral rc * vc
Co   focvxc:(lpert only): integral rc * (dvxc/drho * rho)
Co   rep   :integral rl*exc
Co   repx  :integral rl*ex
Co   repc  :integral rl*ec
Co   rmu   :integral rl*vxc
Co   vl    :local exchange-correlation potential added to vl
Co   fl    :fl(:,:,1:nsp) = vxc by L
Co         :fl(:,:,1+nsp) = exc by L
Co         :fl is only made if 1s digit isw is nonzero
Cl Local variables
Cl   rp    :list of points on the spherical mesh
Cl   exc   :exchange density on (nr,np) mesh
Cl   vxc   :exchange potential on (nr,np) mesh
Cr Remarks
Cr   For perturbation treatment, take numerical derivatives
Cr   df/dr = d/dr (vxc*r**alfa) instead of d/dr (vxc) because
Cr   f is nearly linear for alpha=2/3.
Cr
Cr   In the spin polarized case, the perturbation density rc is not
Cr   spin polarized.  Thus to calc. vxc(rho+drho, m+dm) - vxc(rho,m)
Cr   we use dm=0 and thus drho1 = drho2 = drho/2; thus
Cr     dvxc = lim_drho->0  vxc(rho+drho,rho1+drho/2) - vxc(rho,rho1)
Cu Updates
Cu   01 Dec 13 Enable calls to libxc library.  LDA only for now.
Cu   20 Nov 09 New 20s digit for isw
Cu   02 Apr 09 New option (10's digit isw)
Cu   14 Jun 02 repx and repc (T. Miyake)
Cu    8 Feb 02 focex and focec (T. Miyake)
Cu   13 Jun 00 spin polarized
Cu    3 Jun 00 adapted from old FP code
Cu    10.07.96 dln: modified for Perdew GGA
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nlm,nsp,lpert,lxcf,np,isw
      double precision ri(*),yl(np,nlm),wp(np),rep(2),repx(2),repc(2),
     .  rmu(2),rwgt(nr),focexc(2),focex(2),focec(2),focvxc(2),
     .  rl(nr,nlm,nsp),vl(nr,nlm,nsp),fl(nr,nlm,nsp+1),rc(nr),bscal
      character buf*(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: rp(:,:),rps(:,:,:),exc(:,:),
     .  ex(:,:),ec(:,:),drps(:,:,:)
      real(8), allocatable :: vxc(:,:,:),vx(:,:,:),
     .  vc(:,:,:),vxc2(:,:,:),vx2(:,:,:),vc2(:,:,:)
C ... Local parameters
      integer ilm,ip,ipr,ir,i,n1,ineg(2),isw1,stdo
      double precision weight,suml(0:41),fac,alfa,sum(2)
      double precision rhot,f1,f2,dfdr,dvdr,f
      double precision qs(2),tol,rpneg,xx(1)
      procedure(integer) :: nglob
      procedure(real(8)) :: dmach

      allocate (rp(nr,np),rps(nr,np,nsp),exc(nr,np),ex(nr,np),ec(nr,np))
      allocate (vxc(nr,np,nsp),vx(nr,np,nsp),vc(nr,np,nsp))
      if (lpert /= 0) then
        allocate(vxc2(nr,np,nsp),vx2(nr,np,nsp),vc2(nr,np,nsp))
      endif

      call getpr(ipr)
      stdo = nglob('stdo')

C     'small' density
      tol = 1d-15
      isw1 = mod(isw/10,10)

C     for numerical differentiation
      fac = dmach(1)**(1d0/3d0)
      alfa = 2d0/3d0
      n1 = nr*np
C     call prrmsh('rl',ri,rl,nr,nr,nlm*2)

C --- Generate density point-wise through sphere ---
      call dpzero(rp,n1)
      ineg(1) = 0
      ineg(2) = 0
      rpneg = 0
      rep=0; repx=0; repc=0; rmu=0; qs=0
C     call fyl2p(nr,nlm,nsp,np,yl,rl,rps)
      do  i = 1, nsp
        call dgemm('N','T',nr,np,nlm,1d0,rl(1,1,i),nr,yl,np,0d0,
     .    rps(1,1,i),nr)

C   ... Count and reset to positive points with negative density
        if (isw1 == 1 .or. isw1 == 3) then
          do  ip = 1, np
          do  ir = 1, nr
            if (rps(ir,ip,i) < 0d0) then
              rpneg = min(rpneg,rps(ir,ip,i))
              weight = (ri(ir)**2*rwgt(ir))*wp(ip)
              qs(i)  = qs(i) + rps(ir,ip,i)*weight
              ineg(i) = ineg(i) + 1
              rps(ir,ip,i) = tol
            endif
          enddo
          enddo
        endif
        call daxpy(n1,1d0,rps(1,1,i),1,rp,1)
      enddo

      if (ineg(1)+ineg(2) /= 0 .and. ipr>=10) then
        call awrit5(' vxcnsl (warning) negative rho,'
     .    //' %?#n==2#(%i,%i)#%i%j# points:'//
     .    '  min rho=%;3g   dQ=%;3g',buf,100,stdo,
     .    nsp,ineg(1),ineg(2),rpneg,-qs(1)-qs(2))
      endif

C --- Perturbation treatment: df/dr in vxc2 ---
      if (lpert /= 0 .and. lxcf /= 0) then

C   ... Stop if not CA or BH functional for now
        if (lxcf > 2)
     .    call rx('vxcnsp: Perturbation treatment'//
     .    ' is not implemented for for this functional')

C       Add rp*fac/2 into rp+, rp- and fac*rp into rp
        if (nsp == 2) then
          do  i = 1, nsp
            call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
          enddo
        endif
        call dscal(n1,1+fac,rp,1)
C       Exchange potential at rp+drho
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,1,lxcf,
     .      exc,ex,ec,vxc2(1,1,i),vx2(1,1,i),vc2(1,1,i))
        enddo
C       Restore rp,rps; also add -drho to rp and -drho/2 to rps
        call dscal(n1,1/(1+fac),rp,1)
        if (nsp == 2) then
          do  i = 1, nsp
            call daxpy(n1,-fac,rp,1,rps(1,1,i),1)
          enddo
        endif
        call dscal(n1,(1-fac),rp,1)
C       Exchange potential at rp-drho
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,1,lxcf,
     .      exc,ex,ec,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
        enddo
C       Restore rp,rps
        call dscal(n1,1/(1-fac),rp,1)
        if (nsp == 2) then
          do  i = 1, nsp
            call daxpy(n1,fac/2,rp,1,rps(1,1,i),1)
          enddo
        endif

        do  i = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              rhot = rp(ir,ip)
              if (rhot > 0) then
                f1 = vxc (ir,ip,i)*(rhot*(1-fac))**alfa
                f2 = vxc2(ir,ip,i)*(rhot*(1+fac))**alfa
                dfdr = (f2-f1)/(2d0*fac*rhot)
                vxc2(ir,ip,i) = dfdr
              else
                vxc2(ir,ip,i) = 0
              endif
            enddo
          enddo
        enddo
      else
        call dpzero(ex,n1)
        call dpzero(ec,n1)
      endif  ! End of perturbation treatment

C --- vxc, exc for unperturbed density ---
      if (lxcf > 1000) then ! A libxc functional (must be a local one)
        call xc_libxc(n1,nsp,-lxcf,rps,rps(1,1,nsp),
     .    xx,xx,xx,xx,xx,xx,xx,xx,xx,xx,ex,ec,vx,vc,
     .    vx(1,1,nsp),vc(1,1,nsp),exc,vxc,vxc(1,1,nsp),xx)
C        print *, sngl(exc(100,1)),sngl(vxc(100,1,1)),sngl(vxc(100,1,2))

      elseif (lxcf > 2) then
        i = 1 ; if (nsp == 2) i = 2 ! Enables vxcnsp to pass bounds check
        call evxcp(rps,rps(1,1,nsp),n1,nsp,lxcf,ex,ec,exc,
     .    vx(1,1,1),vx(1,1,i),vc(1,1,1),vc(1,1,i),vxc,vxc(1,1,nsp))

      elseif (lxcf > 0) then
        do  i = 1, nsp
          call evxcv(rp,rps(1,1,i),n1,nsp,lxcf,
     .      exc,ex,ec,vxc(1,1,i),vx(1,1,i),vc(1,1,i))
        enddo
C       print *, sngl(exc(100,1)),sngl(vxc(100,1,1)),sngl(vxc(100,1,2))

      else ! lxcf = 0
        call dpzero(vx,nr*np*nsp)
        call dpzero(vc,nr*np*nsp)
        call dpzero(vxc,nr*np*nsp)
        call dpzero(ex,nr*np)
        call dpzero(ec,nr*np)
        call dpzero(exc,nr*np)
      endif

      if (isw1 == 2) then
        do  i = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              if (rp(ir,ip) <= 0 .or. rps(ir,ip,i) <= 0) then
                vxc(ir,ip,i) = 0
                vx(ir,ip,i) = 0
                vc(ir,ip,i) = 0
              endif
            enddo
          enddo
        enddo
      endif

C --- Optionally scale bxc ---
      call bxcscale(n1,nsp,bscal,n1,vxc)

C --- If rps changed, add Yl-projection of drps into rl ---
      if (ineg(1)+ineg(2) /= 0 .and. isw1 == 3) then
        allocate (drps(nr,np,nsp))
        do  i = 1, nsp
          call dgemm('N','T',nr,np,nlm,-1d0,rl(1,1,i),nr,yl,np,0d0,
     .      drps(1,1,i),nr)
        enddo
        call daxpy(n1*nsp,1d0,rps,1,drps,1)
C       call prrmsh('drps',ri,drps,nr,nr,np*nsp)
C       call prrmsh('rl before',ri,rl,nr,nr,nlm*nsp)
        do   i = 1, nsp
          call dgemm('N','N',nr,nlm,np,1d0,drps(1,1,i),nr,yl,np,1d0,
     .      rl(1,1,i),nr)
        enddo
C       call prrmsh('rl after',ri,rl,nr,nr,nlm*nsp)
        deallocate (drps)
      endif

C --- Scale yl by wp for fast multiplication --- ---
      do  ilm = 1, nlm
        do  ip = 1, np
          yl(ip,ilm) = yl(ip,ilm)*wp(ip)
        enddo
      enddo

      call vxcns4(0,ri,nr,rwgt,np,wp,nsp,rps,exc,ex,ec,vxc,
     .  rep,repx,repc,rmu,qs,buf)

C --- Add perturbation to vxc ---
C     Integrals focexc = int rc vxc, focvxc= rc * dvxc/dr * rhot
      if (lpert /= 0) then
C        if (bscal /= 1d0 .and. nsp == 2) then
C          call rx('pert form of vxc not implemented with bscal')
C        endif
        focvxc(1) = 0
        focvxc(2) = 0
        focexc(1) = 0
        focexc(2) = 0
        focex(1)  = 0
        focex(2)  = 0
        focec(1)  = 0
        focec(2)  = 0
        do  i  = 1, nsp
          do  ip = 1, np
            do  ir = 1, nr
              rhot = rp(ir,ip)
              if (rhot <= 0 .or. rps(ir,ip,i) <= 0) then
                vxc2(ir,ip,i) = 0
              endif
C             Debugging
C              if (rps(ir,ip,i) < 0 .or. rhot < 0) then
C                if (vxc(ir,ip,i) /= 0 .or. vxc2(ir,ip,i) /= 0) then
C                  print *, vxc(ir,ip,i), vxc2(ir,ip,i)
C                  stop 'oops'
C                endif
C              endif
              if (rps(ir,ip,i) > 0 .and. rhot > 0) then
                f  = vxc(ir,ip,i)*rhot**alfa
                dfdr = vxc2(ir,ip,i)
                dvdr = (dfdr - alfa*f/rhot) / rhot**alfa
                weight = (ri(ir)**2*rwgt(ir))*wp(ip) * rc(ir)
                focvxc(i) = focvxc(i) + weight*dvdr*rps(ir,ip,i)
                focexc(i) = focexc(i) + weight*vxc(ir,ip,i)/nsp
                focex(i)  = focex(i)  + weight*vx(ir,ip,i)/nsp
                focec(i)  = focec(i)  + weight*vc(ir,ip,i)/nsp
                vxc(ir,ip,i) = vxc(ir,ip,i) + dvdr*rc(ir)
              endif
            enddo
          enddo
        enddo
      endif

      if (lxcf == 0) goto 99

C --- Add Yl-projection of vxc into vl ---
      do  i = 1, nsp
        call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,1d0,
     .    vl(1,1,i),nr)
      enddo

C --- Optionally make l-resolved vxc_L,exc_L ---
      if (mod(isw,10) == 1) then
        do  i = 1, nsp
          call dgemm('N','N',nr,nlm,np,1d0,vxc(1,1,i),nr,yl,np,0d0,
     .      fl(1,1,i),nr)
        enddo
        call dgemm('N','N',nr,nlm,np,1d0,exc,nr,yl,np,0d0,
     .    fl(1,1,1+nsp),nr)
      endif

C --- Print out int (rl*vl) resolved by l ---
      if (ipr > 30) then
        call vxcns5(0,ipr,'rho*vtot',nlm,nsp,nr,ri,rwgt,rl,vl,suml,sum,buf)
C      lmax = ll(nlm)
C      do  42  i = 1, nsp
C      do  43  l = 0, lmax
C   43 suml(l) = 0d0
C      do  40  ilm = 1, nlm
C      l = ll(ilm)
C      do  40  ir = 1, nr
C   40 suml(l) = suml(l) + rl(ir,ilm,i)*vl(ir,ilm,i)*ri(ir)**2*rwgt(ir)
C      if (i == 1) write(stdo,341) (suml(l),l = 0,lmax)
C      if (i == 2) write(stdo,342) (suml(l),l = 0,lmax)
C  341 format(' rho*vxc by l: ',f13.6,4f10.6:/(18x,4f10.6))
C  342 format('       spin 2: ',f13.6,4f10.6:/(18x,4f10.6))
C   42 continue
      endif

   99 continue
      deallocate (rp,rps,exc,ex,ec)
      deallocate (vxc,vx,vc)
      if (lpert /= 0) then
        deallocate(vxc2,vx2,vc2)
      endif

      end
      subroutine vxcns3(nr,nlm,nsp,ri,rl,isgn)
C- Scales rho by r**2, or undoes scaling
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlm   :L-cutoff for density expansion
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ri    :radial mesh points
Ci   isgn  :1, scale rl by 1/r**2
Ci         :else scale rl by r**2
Co Outputs
Co   rl   :rl is scaled by r**2 or 1/r**2, depending on isgn
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer nr,nlm,nsp,isgn,i,ilm,ir
      double precision rl(nr,nlm,nsp),ri(nr),rho2,rho3

C --- Scale rho by 1/r**2 ---
      if (isgn == 1) then
        do  i = 1, nsp
          do  ilm = 1, nlm
            rl(1,ilm,i) = 0d0
            do  ir = 2, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)/ri(ir)**2
            enddo
          enddo
        enddo
C  ...  Extrapolate rho to origin
        do  i = 1, nsp
        rho2 = rl(2,1,i)
        rho3 = rl(3,1,i)
          rl(1,1,i) = (rho2*ri(3)-rho3*ri(2))/(ri(3)-ri(2))
        enddo
      else
        do  i = 1, nsp
          do  ilm = 1, nlm
            do  ir = 1, nr
              rl(ir,ilm,i) = rl(ir,ilm,i)*ri(ir)**2
            enddo
          enddo
        enddo
      endif
      end

      subroutine vxcns4(isw,ri,nr,rwgt,np,wp,nsp,rp,exc,ex,ec,vxc,
     .  rep,repx,repc,rmu,qs,buf)
C- Integrals of vxc, reps and rmu in sphere.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcf  :if not positive, return qs only
Ci   isw   :10s digit
Ci         : 2 for any point where rho<0 or rho_isp<0, zero potential
Ci   ri    :radial mesh points
Ci   nr    :number of radial mesh points
Ci   rwgt  :radial mesh weights
Ci   np    :number of spherical mesh points
Ci   wp    :spherical mesh weights
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cl   rp    :list of points on the spherical mesh
Cl   ex    :exchange             energy density on (nr,np) mesh
Cl   ec    :correlation          energy density on (nr,np) mesh
Cl   exc   :exchange-correlation energy density on (nr,np) mesh
Cl   vxc   :exchange correlation potential on (nr,np) mesh
Co Outputs
Co   qs    :spin-resolved charge
Co   ... The following are returned if lxcf>0
Co   rep   :integral rl*exc
Co   repx  :integral rl*ex
Co   repc  :integral rl*ec
Co   rmu   :integral rl*vxc
Cl Local variables
Cl   rpneg :number of points at which density < 0
Cr Remarks
Cu Updates
Cu   12 Mar 04 created from vxcns2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,np,isw
      double precision ri(*),wp(np),rep(2),repx(2),repc(2),rmu(2),
     .  rwgt(nr),rp(nr,np,nsp),exc(nr,np),ex(nr,np),ec(nr,np),vxc(nr,np,nsp)
      character*(*) buf
C ... Local parameters
      integer ip,ipr,ir,i,stdo,isw1
      double precision weight,rpneg,qs(2)
      procedure(integer) :: lgunit

      stdo = lgunit(1)
      call getpr(ipr)
      isw1 = mod(isw/10,10)

C     call awrit0('%&',' ',0,0)

C --- Integrals reps, rmu ---
      rpneg = 0
      do  i = 1, nsp
        qs(i) = 0d0
        rep(i) = 0d0
        repx(i) = 0d0
        repc(i) = 0d0
        rmu(i) = 0d0
        do  ip = 1, np
        do  ir = 1, nr
          rpneg = min(rpneg,rp(ir,ip,i))
          weight = (ri(ir)**2*rwgt(ir))*wp(ip)
          qs(i)  = qs(i)  + rp(ir,ip,i)*weight
C         if (lxcf == 0) cycle
          if (isw1 == 2) then
          if (rp(ir,ip,1)+rp(ir,ip,nsp) <= 0 .or. rp(ir,ip,i) < 0) then
            exc(ir,ip) = 0
            ex(ir,ip) = 0
            ec(ir,ip) = 0
            vxc(ir,ip,i) = 0
          endif
          endif
C         Debugging
C          if (rp(ir,ip,1) < 0 .or. rp(ir,ip,nsp) < 0) then
C            if (vxc(ir,ip,i) /= 0) then
C              print *, vxc(ir,ip,i)
C              stop 'oops'
C            endif
C          endif
          rep(i) = rep(i) + exc(ir,ip)*rp(ir,ip,i)*weight
          repx(i) = repx(i) + ex(ir,ip)*rp(ir,ip,i)*weight
          repc(i) = repc(i) + ec(ir,ip)*rp(ir,ip,i)*weight
          rmu(i) = rmu(i) + vxc(ir,ip,i)*rp(ir,ip,i)*weight
        enddo
        enddo
C        if (ipr >= 30 .and. i == 1) write(stdo,725) rmu(i),rep(i),qs(i)
C        if (ipr >= 30 .and. i == 2) write(stdo,726) rmu(i),rep(i),qs(i),
C     .    rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2)
C  725   format(' vxcnsp: loc rmu=',f11.6,'  rep=',f11.6,'  q = ',f10.6)
C  726   format(' spin 2:         ',f11.6,'      ',f11.6,'      ',f10.6/
C     .         '  total:         ',f11.6,'      ',f11.6,'      ',f10.6)

        if (ipr >= 30 .and. i == 1) call awrit3(' vxcnsp: loc rmu=%;11,6D  rep=%;11,6D  q =%;11,6D',
     .                              buf,100,stdo,rmu(i),rep(i),qs(i))
        if (ipr >= 30 .and. i == 2) call awrit6(' spin 2:         %;11,6D      %;11,6D     %;11,6D%N'//
     .                                          '  total:         %;11,6D      %;11,6D     %;11,6D',
     .    buf,140,stdo,rmu(i),rep(i),qs(i),rmu(1)+rmu(2),rep(1)+rep(2),qs(1)+qs(2))

      enddo

      if (rpneg < 0 .and. ipr >= 10) call awrit1(
     .  ' vxcnsp (warning): negative rho: min val = %;3g',buf,100,stdo,rpneg)

      end

      subroutine vxcns5(isw,ipr,strn,nlml,nsp,nr,ri,rwgt,rl,vl,suml,sum,buf)
C- Integrals of rl*vl resolved by l
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   :0 rl is true density
Ci         :1 rl is true density * ri**2
Ci   strn  :used in printout
Ci   ri    :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :number of Ylm's in rl and vl
Co Outputs
Co   suml  :l-and spin-resolved integrals of rl*vl
Co   sum   :spin-resolved integrals of rl*vl
Cr Remarks
Cu Updates
Cu   12 Mar 04 created from vxcns2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character strn*(*),buf*(*)
      integer isw,ipr,nlml,nsp,nr
      double precision ri(nr),rwgt(nr),rl(nr,nlml,nsp),vl(nr,nlml,nsp),
     .  suml(0:20,2),sum(2)
C ... Local parameters
      integer ir,stdo,lgunit,lmax,ilm,isp,l,ll
      double precision dsum,dot3

C     call awrit0('%&',' ',0,0)

      stdo = lgunit(1)
      lmax = ll(nlml)

      do   isp = 1, nsp
      call dpzero(suml(0,isp),lmax+1)
      do  ilm = 1, nlml
        l = ll(ilm)
        if (isw == 1) then
        suml(l,isp) = suml(l,isp) +
     .    dot3(nr,rl(1,ilm,isp),vl(1,ilm,isp),rwgt)
        else
        do  ir = 1, nr
        suml(l,isp) = suml(l,isp) +
     .                rl(ir,ilm,isp)*vl(ir,ilm,isp)*ri(ir)**2*rwgt(ir)
        enddo
        endif
      enddo

      sum(isp) = dsum(lmax+1,suml(0,isp),1)

      if (ipr > 30) then
C        if (isp == 1) write(stdo,341) strn,sum(isp),(suml(l,isp),l = 0,lmax)
C        if (isp == 2) write(stdo,342) sum(isp),(suml(l,isp),l = 0,lmax)

        if (isp == 1) call awrit4(' '//strn(1:8)//'%;13,6D ... by l: %;12,6D%n;10,6D',
     .    buf,100,stdo,sum(isp),suml(0,isp),lmax,suml(1,isp))
        if (isp == 2) call awrit4('  spin 2:%;13,6D ... by l: %;12,6D%n;10,6D',
     .    buf,100,stdo,sum(isp),suml(0,isp),lmax,suml(1,isp))

  341   format(1x,a8,f13.6,' ... by l: ',f12.6,4f10.6:/(18x,4f10.6))
  342   format('  spin 2:',f13.6,' ... by l: ',f12.6,
     .    4f10.6:/(18x,4f10.6))
      endif

      enddo

      end
