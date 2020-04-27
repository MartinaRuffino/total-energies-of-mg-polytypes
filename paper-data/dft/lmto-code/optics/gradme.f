      subroutine gradme(opt,nf1,nf2,isp,nsp,nspc,ic1,ic2,lmx,lmxa,nl2,rgrad,gradm,zgradm)
C- Matrix elements <phi grad phi'> from radial matrix elements
C ----------------------------------------------------------------
Ci Inputs
Ci   opt   :passed as mode to grmes, which see
Ci         :  0 return gradients grad(x,y,z) in real harmonic basis
Ci         :  1 return gradients grad(x,y,z) in spherical harmonic basis, m=-l...l order
Ci         :  2 return gradients grad(x,y,z) in spherical harmonic basis, m=l...-l order
Ci         :  3 return gradients grad(x,y,z) in spherical harmonic basis, old conventions
Ci         :    DOES NOT WORK.
Ci         :10s digit controls form of gradm returned by grmes
Ci         :  0 Return real part in gradm, discarding the imaginary part
Ci         :    Often (e.g. for opt=0 or opt=201 or 202) the result should be real anyway.
Ci         :    In this case gradm = gradm(k,n2,n1,ilm2,ilm1), as discussed in Remarks
Ci         :  1 Return gradients as COMPLEX matrix zgradm, with indices same as opt/10=0.
Ci         :  2 Same as opt/10=0, but indices are permuted; see pgrmes
Ci         :  3 Same as opt/10=1, but indices are permuted; see pgrmes
Ci         :100s digit type of gradient
Ci         :  0 return gradients grad(x,y,z)
Ci         :  1 return gradients grad(+,-,z)
Ci         :  Add 2 to swap components 1 and 2 (x <-> y or + <-> -)
Ci         :Use 101 for SH basis, (+,-,z) gradient and 301 for SH basis, (-,+,z) gradient
Ci OLD
Ci         :1s digit opt controls whether spherical or real harmonics
Ci         :You can choose gradm returned in real or spherical coordinates,
Ci         :and also whether grad_m corresponds to Cartesian or spherical vectors.
Ci         :0 return gradients grad(x,y,z) in real harmonic basis
Ci         :1 return gradients grad(x,y,z) in spherical harmonic basis
Ci         :2 rotate grad(x,y,z) to grad(+,-,z)
Ci         :4 swap components 1 and 2 (x <-> y or + <-> -)
Ci         :  Use 3 for SH basis, (+,-,z) ordering and 7 for SH basis, (-,+,z) ordering.
Ci
Ci         :10s digit opt = even => return REAL gradm; odd => return COMPLEX zgradm
Ci         :0 Return real part in gradm, discarding the imaginary part
Ci         :  For "normal" modes (opt=0, 3 or 7) the result should be real anyway.
Ci         :  In this case gradm = gradm(k,n2,n1,ilm2,ilm1), as discussed in Remarks
Ci         :1 Return gradients as COMPLEX matrix zgradm.
Ci         :2 Same as opt/10=0, but indices are permuted; see gradm below
Ci         :3 Same as opt/10=1, but indices are permuted; see gradm below
Ci   nf1   :Number of gradient functions
Ci   nf2   :Number of pair functions
Ci   isp   :current spin channel for collinear spin polarized case
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ic1,ic2 range of classes or sites for which to make matrix elements
Ci   lmx   :dimensions rgrad integrals from 0 to lmx
Ci   lmxa  :make integrals for each ic=(ic1:ic2) from 0 to lmxa(ic)
Ci   nl2   :dimensions gradm
Ci   rgrad(if1,if2,ll,l,nsp):  <g2_l' grad_r g1_l> with l' = l +/- 1
Ci         :if1 = index to ket radial partial wave g1
Ci         :if2 = index to bra radial partial wave g2
Ci         :Example: if1 = 1,2 for phi,phidot:
Ci         :    if1   if2  g1      g2
Ci         :     1     1   phi     phi
Ci         :     2     1   phidot  phi
Ci         :     1     2   phi     phidot
Ci         :     2     2   phidot  phidot
Ci         :ll= 1 => l' = l+1  <g2(l') | grad_r g1(l)> - (l+1) < g2(l') | 1/r g1(l) >
Ci         :    2 => l' = l-1  <g2(l') | grad_r g1(l)> +    l  < g2(l') | 1/r g1(l) > 
Co Outputs
Co   gradm(k,if1,if2,ilm2,ilm1,isp,ic):
Co         :<g2 grad g1> for g1 = g1(r) Y_ilm1(rhat) and g2 = g2(r) Y_ilm2(rhat)
Co         :1s bit opt 0:  Ylm are in real harmonics
Co         :1s bit opt 1:  Ylm are in spherical harmonics
Co         :2s bit opt 0:  k=1 : grad_x  k=2 : grad_y   k=3 : grad_z
Co         :2s bit opt 1:  k=1 : grad_+  k=2 : grad_-   k=3 : grad_z
Co         :If 4s bit opt is set, swap k=1 with k=2.
Co         :If 10s digit opt>2 permute gradm(k,if1,if2,ilm2,ilm1) -> gradm(,ilm2,ilm1,k,if1,if2)
Co   zgradm:complex version of gradm
Co         :Whether gradm or zgradm is returned depends on 10s digit opt
Co   Noncollinear case: gradm = gradm(k,if1,if2,ilm1,ilm2,isp,jsp,ic)
Cl Local variables
Cl   ksp   : isp  in the collinear case, isp = 1 ... nsp
Cl         : ispc in the noncollinear case
Cl   kspc  : isp  in the collinear case, isp = 1 ... nsp
Cl         :  compound spin index (isp,ispc) in the noncollinear case
Cr Remarks
Cr   Adapted from optics package originally developed by V. Antropov
Cr   User chooses representation of gradient operator; see opt above and also
Cr   questaal.org/docs/numerics/spherical_harmonics/#gradients-of-spherical-harmonics
Cu Updates
Cu   05 May 18 Code rewritten from scratch.
Cu             Sign convention for grad(k=2) is flipped from prior version
Cu   22 Jun 14 Extended to the noncollinear case
Cu   30 Dec 13 Redesign to use with FP package
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,isp,nsp,nspc,nf1,nf2,nl2,ic1,ic2,lmx,lmxa(ic2)
      double precision rgrad(nf1,nf2,2,0:lmx,nsp*nspc,ic2)
      double precision gradm(3,nf1,nf2,nl2,nl2,nsp*nspc,ic2)
      double complex  zgradm(3,nf1,nf2,nl2,nl2,nsp*nspc,ic2)
C ... Local parameters
      integer ic,ispc,ksp,kspc
C#ifdefC DEBUG
C      integer il,ilmun,lun,mun,ilmoc,loc,moc,if1,if2,lopt
C      double precision pp,pm,pz
CC     double precision grads(3,nf1,nf2,nl2,nl2),err
C      double precision xgradm(nl2,nl2,3,nf1,nf2,nsp*nspc,ic2)
C      procedure(real(8)) :: dlength
C
C      goto 100
C#endif

C --- Loop over classes or sites, and noncollinear spins ---
      do  ispc = 1, nspc
        ksp = isp ; if (nspc == 2) ksp = ispc
        kspc = isp + 2*(ispc-1)

        do  ic = ic1, ic2
          call grmes(opt,lmxa(ic),rgrad(1,1,1,0,kspc,ic),nf1,nf2,[0],[0],nl2,
     .      gradm(1,1,1,1,1,kspc,ic),zgradm(1,1,1,1,1,kspc,ic))
C         call yprmi('gradm ic=%i opt=%i',ic,opt,1,gradm(1,1,1,1,1,ic,1),0,3*nf1*nf2,3*nf1*nf2,nl2*nl2)
        enddo                   ! Loop over classes or sites
      enddo                     ! Loop over noncollinear spins


C --- Matrix elements <g1 grad g2> using Wigner-Eckhart theorem ---
C#ifdefC DEBUG
C  100 continue
C      call dpzero(xgradm,size(xgradm))
C
C      do  ispc = 1, nspc
C      ksp = isp ; if (nspc == 2) ksp = ispc
C      kspc = isp + 2*(ispc-1)
C      do  ic = ic1, ic2
C
CC       For each ilmun, lmuc pair, do
C        ilmun = 0
C        do  lun = 0, lmxa(ic)
C          do  mun = -lun, lun
C            ilmun = ilmun+1
C            ilmoc = 0
C            do  loc = 0, lmxa(ic)
C              il = 1
C              if (loc-lun == 1) il = 2
C              do  moc = -loc, loc
C                ilmoc = ilmoc+1
C
CC               +,-,z matrix elements for each function pair
C                outer : do  if2 = 1, nf2
C                do  if1 = 1, nf1
C                  call pgrdme(lun,mun,loc,moc,
C     .              rgrad(if1,if2,il,loc,kspc,ic),
C     .              rgrad(if1,if2,il,loc,kspc,ic),pp,pm,pz)
C                  xgradm(ilmun,ilmoc,1,if2,if1,kspc,ic) = pp
C                  xgradm(ilmun,ilmoc,2,if2,if1,kspc,ic) = pm
C                  xgradm(ilmun,ilmoc,3,if2,if1,kspc,ic) = pz
C
C
C                  exit outer
C
C                enddo
C                enddo outer
C
C              enddo
C            enddo
C          enddo
C        enddo
C
C        if1 = 1 ; if2 = 1
C        call prmx('ycg+',xgradm(1:nl2,1:nl2,1,if1,if2,1,ic),nl2,nl2,nl2)
C        call prmx('ycg-',xgradm(1:nl2,1:nl2,2,if1,if2,1,ic),nl2,nl2,nl2)
C        stop
C
C      enddo  ! Loop over classes or sites
C      enddo  ! Loop over noncollinear spins
C
C#endif

      end

      subroutine pgrdme(lp,mp,l,m,rp,rm,pp,pm,pz)
C- Kernel called by gradme
C ----------------------------------------------------------------
Ci Inputs
Ci   lp    :l' which should be l+1 or l-1
Ci   mp    :m' which should be m+1 or m-1 or m
Ci   l,m   :index to Y_lm
Ci   rp    :Radial matrix element (l+1 | grad_r | l)
Ci   rm    :Radial matrix element (l-1 | grad_r | l)
Co Outputs
Co   Refer to coefficients A on the Questaal web page:
Co   questaal.org/docs/numerics/spherical_harmonics/#gradients-of-spherical-harmonics
Co                          coff   case       coff   case       else
Co   pp    : grad- mp=m+1   A+_+1  lp=l+1     A-_+1  lp=l-1     0
Co   pm    : grad+ mp=m-1   A+_-1  lp=l+1     A-_-1  lp=l-1     0
Co   pz    : gradz mp=m     A+_0   lp=l+1     A-_0   lp=l-1     0
Co   See Remarks
Cr Remarks
Cr   Adapted from Edmonds, "Angular Momentum in Quantum Mechanics" with sign errors corrected
Cr   Sign checks for A+ against Edmonds Section 5.7
Cr   Each r^l*Ylm is of the form (polynomial in (rp,rm,z) * radial factor
Cr   Assume r^l*Ylm are even functions of rp = (-x-iy) and rm = (x-iy) (true by inspection for l=0,1,2,3)
Cr
Cr   Consider raising operator: use phi(r) = r^-l-1  => dphi/dr = (-l-1) phi/r
Cr   Note that d(1/r^n)/d(rp,rm,z) is (positive,positive,negative)
Cr   1. Erroneous minus sign for A+0.
Cr      r^1 grad_z Y00/r   = -1/sqrt(3)   Y(1,0)
Cr      r^3 grad_z Y10/r^2 = -2*sqrt(3/5) Y(2,0)
Cr      r^5 grad_z Y20/r^3 = -3*sqrt(5/7) Y(3,0)
Cr   All are negative because d(1/r^n)/dz < 0.
Cr   Sign of A+0 should always be positive since radial factor is negative
Cr   No factor (-1)^(l+m)
Cr   2. Check A+_(-1,+1)
Cr      Derivatives wrt rp, rm of product  1/r^n Y_lm yield two both positive.
Cr      Factor dphi/dr is negative.
Cr   Sign of A+_(-1,1) should always be negative since radial factor is negative
Cr   No factor (-1)^(l+m)
Cr   Changes 1 and 2 yield exact agreement with grmes.
Cr
Cr   Lowering operator: Consider phi(r) = r^l  => dphi/dr + (l+1)phi/r = (2*l+1) phi/r
Cr   1. Check (-1)^(l+m) for Y_l0.  By inspection, grad_z Y_l0 is always positive.
Cr      Factor (-1)^(l+m) should not be present for A-_0.
Cr   2. Y_1m, Y2m, Y3m are even functions of rp = (-x-iy) and rm = (x-iy) by inspection
Cr      Derivatives should contain no (-) signs.
Cr      Sign of A-_(-1,1) should always be positive since radial factor is positive
Cr
Cr   These signs are consistent with Varshalovich "Quantum Theory of Angular Momentum" 5.8.3,
Cr   except for factor of -1 for (+,-) components?
Cr
Cr   Another apparent inconsistency with Edmonds and with Varshalovich:
Cr   Consider Y1(-1,1)*r = sqrt(3)*Y00(rm,rp)
Cr     grad(-1) Y1(-1)*r = sqrt(3)*Y00 : l->l-1 and m->m+1.  Corresponds to A-_+1
Cr     grad(+1) Y1(+1)*r = sqrt(3)*Y00 : l->l-1 and m->m-1.  Corresponds to A-_-1
Cr   Consider Y00/r
Cr     grad(-1) Y00/r propto rp propto Y11  : l->l+1 and m->m+1.  Corresponds to A+_+1
Cr     grad(+1) Y00/r propto rm propto Y1-1 : l->l+1 and m->m-1.  Corresponds to A+_-1
Cr   Thus what both Edmonds and Varshalovich call grad_(+1,-1) are actually grad_(-1,+1)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lp,mp,l,m
      double precision rp,rm,pp,pm,pz
C ... Local parameters
      double precision r,d

      pp = 0d0; pm = 0d0; pz = 0d0
      if (abs(m) > l .or. abs(lp-l) /= 1) return

      if (lp == l-1) then
        r=rm; d = 2*(2*l+1)*(2*l-1)
      else
        r=rp; d = 2*(2*l+3)*(2*l+1)
      endif

C      As written in Edmonds Section 5.7
C      m1 = (-1)**(l+m)
CC     m1 = 1
C
C      if (lp == l+1 .and. mp == m-1) pm = m1*dsqrt((l-m+1)*(l-m+2)/d)*r    ! A+(-1)
C      if (lp == l+1 .and. mp == m+1) pp = m1*dsqrt((l+m+1)*(l+m+2)/d)*r    ! A+(1)
C      if (lp == l+1 .and. mp == m)   pz = -m1*dsqrt(2*(l+m+1)*(l-m+1)/d)*r ! A+(0)
C
C      if (lp == l-1 .and. mp == m-1) pm = m1*dsqrt((l+m-1)*(l+m)/(d))*r    ! A-(-1)
C      if (lp == l-1 .and. mp == m+1) pp = m1*dsqrt((l-m-1)*(l-m)/(d))*r    ! A-(1)
C      if (lp == l-1 .and. mp == m)   pz = m1*dsqrt(2*(l+m)*(l-m)/d)*r      ! A-(0)

      if (lp == l+1 .and. mp == m-1) pm =-dsqrt((l-m+1)*(l-m+2)/d)*r    ! A+(-1)
      if (lp == l+1 .and. mp == m+1) pp =-dsqrt((l+m+1)*(l+m+2)/d)*r    ! A+(1)
      if (lp == l+1 .and. mp == m)   pz = dsqrt(2*(l+m+1)*(l-m+1)/d)*r  ! A+(0)

      if (lp == l-1 .and. mp == m-1) pm = dsqrt((l+m-1)*(l+m)/(d))*r    ! A-(-1)
      if (lp == l-1 .and. mp == m+1) pp = dsqrt((l-m-1)*(l-m)/(d))*r    ! A-(1)
      if (lp == l-1 .and. mp == m)   pz = dsqrt(2*(l+m)*(l-m)/d)*r      ! A-(0)

      end

C#ifdefC TEST
C      subroutine opgrdme(lp,mp,l,m,r1,r2,pp,pm,pz)
CC- Original Kernel called by gradme, written by V. Andropov
CC ----------------------------------------------------------------
CCi Inputs
CCi   lp    :l' which should be l+1 or l-1
CCi   mp    :m' which should be m+1 or m-1 or m
CCi   l,m   :index to grad Y_lm
CCi   r1    :Radial matrix element connecting l to l+/-1
CCi   r2    :Radial matrix element connecting l to l+/-1
CCo Outputs
CCo   pp    :
CCo   pm    :
CCo   pz    :
CCr Remarks:
CCr   sign convention for pm (ME of grad+) is opposite from pgrdme
CC ----------------------------------------------------------------------
C      implicit none
C      integer lp,mp,l,m
C      double precision r1,r2
C      double precision pp,pm,pz
C
C      pp = 0d0
C      pm = 0d0
C      pz = 0d0
C      if (abs(m) > l .or. abs(lp-l) /= 1) return
C
CC     if (lp == l+1 .and. mp. eq. m+1) pp=r1
C      if (lp == l+1 .and. mp == m+1)
C     .  pp = -dsqrt((1+l+m)*(2+l+m)/(2d0*(2*l+1)*(2*l+3)))*r1
C
CC     if (lp == l-1 .and. mp == m+1) pp=r2
C      if (lp == l-1 .and. mp == m+1)
C     .  pp = dsqrt((l-m-1)*(l-m)/(2d0*(2*l-1)*(2*l+1)))*r2
C
CC     if (lp == l+1 .and. mp == m-1) pm=r1
C      if (lp == l+1 .and. mp == m-1)
C     .  pm = dsqrt((l-m+1)*(l-m+2)/(2d0*(2*l+1)*(2*l+3)))*r1
C
CC     if (lp == l-1 .and. mp == m-1) pm=r2
C      if (lp == l-1 .and. mp == m-1)
C     .  pm = -dsqrt((l+m-1)*(l+m)/(2d0*(2*l-1)*(2*l+1)))*r2
C
CC     if (lp == l+1 .and. mp == m) pz=r1
C      if (lp == l+1 .and. mp == m)
C     .  pz = dsqrt((l-m+1)*(l+m+1)/((2d0*l+1)*(2*l+3)))*r1
C
CC     if (lp == l-1 .and. mp == m) pz=r2
C      if (lp == l-1 .and. mp == m)
C     .  pz = dsqrt((l-m)*(l+m)/((2d0*l+1)*(2*l-1)))*r2
C
C       end
C#endif

C#ifdefC TEST
C      subroutine fmain
C      implicit none
C
C      character*(120) strn,cpmz*3,dc*1
C      integer, parameter :: nnn=300
C      integer mode,modeg,nlm,nlmf,lmax,j,ilm,l,lav,m,nf1,nf2,
C     .  ilmun,lun,mun,ilmoc,loc,moc,if1,if2,il,mylm,ifi
C      double precision pi,srfpi,err(3),rp,rm,pp,pm,pz
C      real(8), allocatable :: rgrad(:,:),gradm(:,:,:,:,:),xgradm(:,:,:,:,:)
C      integer, parameter :: NULLI=-99999
C      procedure(logical) :: cmdopt
C      procedure(integer) :: ll,fopng,cmdoptsw,a2vec
CC     data cpmz / 'xyz'/
C      data cpmz / '-+z'/
C
CC     mode=4  => decompose r^-l grad*(r^l Ylm) into linear combination of Ylm, Ylm = sph. harm. (no r^l)
CC     mode=14 => decompose      grad*(Ylm) into linear combination of Ylm
CC     mode=24 => decompose r^(l+1) grad*(r^(-l-1) Ylm) into linear combination of Ylm
CC     mode=34 => read radial matrix elements from file
C
C      mode =  4
CC     mode = 14
C      mode = 24
CC      mode = 34
C      if (cmdopt('--mode=',NULLI,0,strn)) then
C        dc = '='
C        j = 7
C        if (a2vec(strn,len_trim(strn),j,2,', '//dc,3,2,1,ilm,mode) < 1) call rx('failed to parse '//trim(strn))
C      endif
C      if (mode/=4 .and. mode/=14 .and. mode/=24 .and. mode/=34) call rx('bad mode')
C
C
C      lmax = 4
CC     rgrme file should contain lmax first line, then rgrme in free format
C      if (mode == 34) then
C        ifi = fopng('rgrme.in',-1,0)
C        read(ifi,*) lmax
C        print *, 'reading radial matrix elements from file rgrme.in, lmax=',lmax
C      endif
C
C      nlmf = (lmax+2)**2
C      nlm  = (lmax+1)**2
C      pi = 4*datan(1d0)
C      srfpi = dsqrt(4*pi)
C
C      modeg = 321               ! modeg = 27
CC      modeg = 121               ! modeg = 23
CC      modeg = 220               ! modeg = 24
CC      modeg = 20                ! modeg = 20
CC      modeg = 322               ! reverse l,m order
C
C      if (modeg/100 == 0) cpmz = 'xyz'
C      if (modeg/100 == 1) cpmz = '+-z'
C      if (modeg/100 == 2) cpmz = 'yxz'
C      if (modeg/100 == 3) cpmz = '-+z'
C      mylm = 1; if (mod(modeg,10) == 2) mylm = -1
C
C      if (mod(modeg,10) > 2) stop 'oops'
C
C      call info2(1,0,0,' grmes-mode=%i '//
C     .  '%-1j%?#(n%10)%4==0#(real harmonics##'//
C     .  '%-1j%?#(n%10)%4==1#(spherical harmonics (-l:l)##'//
C     .  '%-1j%?#(n%10)%4==2#(spherical harmonics (l:-l)##'//
C     .  ', '//cpmz//' grad)  mode=%i'//
C     .  '%-1j%?#n==4# (gradients of r^l Y_L)##'//
C     .  '%-1j%?#n==14# (gradients of Y_L)##'//
C     .  '%-1j%?#n==34# (radial gradients read from file rgrme.in)##'//
C     .  '%-1j%?#n==24# (gradients of r^-l-1 Y_L)##',modeg,mode)
C
CC     Same result with call to grmes
C      nf1 = 1; nf2 = 1
C      allocate(gradm(nlm,nlm,3,nf1,nf2),xgradm(nlm,nlm,3,nf1,nf2))
C
C      allocate(rgrad(2,0:lmax)); call dpzero(rgrad,size(rgrad))
C      do  l = 0, lmax-1
C        if (mode ==  4) rgrad(2,l) =  2*l+1   ! for (l-1, l)
C        if (mode == 14) rgrad(1,l) =  -l      ! for (l+1, l)
C        if (mode == 14) rgrad(2,l) =  l+1     ! for (l-1, l)
C        if (mode == 24) rgrad(1,l) = -2*l-1   ! for (l+1, l)
C      enddo
C      if (mode == 34) then
C        read(ifi,*) rgrad
C        print *, 'read radial matrix elements from file rgrme.in'
C      endif
C
C      call info0(1,1,0," ... matrix elements from grmes."//"%N%19fl  m    coff ...")
C      call grmes(modeg,lmax,rgrad,1,1,[0],[0],nlm,gradm)
CC     call grmes(0,lmax,rgrad,1,1,[0],[0],nlm,gradm)
C
C      err = 0
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = mylm*(ilm-lav)
C        call awrit2(' '//cpmz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        call pryl(0,trim(strn),ilm,1-mylm,gradm(1:lmax**2,ilm,j,1,1),lmax**2,1d0)
C      enddo
C      print *
C      enddo
C
C      if (modeg /= 321) call cexit(1,1)
CC     if (modeg /= 321 .and. modeg /= 322) call cexit(1,1)
C
C      call info0(1,1,0," ... Compare DIFFERENCES to matrix elements from pgrdme"//"%N%19fl  m    coff ...")
C
CC     For each ilmun, lmuc pair, do
C      ilmun = 0
C      do  lun = 0, lmax
C        do  mun = -lun, lun
C          ilmun = ilmun+1
C          ilmoc = 0
C          do  loc = 0, lmax
C            il = 1
C            if (loc-lun == 1) il = 2
C            do  moc = -loc, loc
C              ilmoc = ilmoc+1
C
CC             +,-,z matrix elements for each function pair
C              if1=1; if2 = 1
C              rp = rgrad(1,loc)
C              rm = rgrad(2,loc)
C
C              call pgrdme(lun,mun,loc,moc,rp,rm,pp,pm,pz)
CC             call opgrdme(lun,mun,loc,moc,rp,rm,pp,pm,pz)
C              xgradm(ilmun,ilmoc,1,if2,if1) = pp  ! ME of grad-
C              xgradm(ilmun,ilmoc,2,if2,if1) = pm  ! ME of grad+
C              xgradm(ilmun,ilmoc,3,if2,if1) = pz  ! MD pf gradz
C
C            enddo               ! moc
C          enddo                 ! loc
C        enddo                   ! nun
C      enddo                     ! lun
C
C      err = 0
C      do  j = 1, 3
C      do  ilm = 1, 16
C        l = ll(ilm)
C        lav = l*l+l+1
C        m = ilm-lav
C        call awrit2(' '//cpmz(j:j)//' grad Yl(%i%,2i): ',strn,len(strn),0,l,m)
C        call pryl(0,trim(strn),ilm,0,xgradm(1:lmax**2,ilm,j,1,1)-1*gradm(1:lmax**2,ilm,j,1,1),lmax**2,1d0)
C        err(j) = err(j) + sum(abs(xgradm(1:lmax**2,ilm,j,1,1)-gradm(1:lmax**2,ilm,j,1,1)))
C      enddo
C      print *
C      enddo
C      call info2(1,0,1,' cumulative diff between (grme, pgrdme) ME ('//cpmz//') = %3:1;3,3g',err,2)
C
C      call cexit(1,1)
C      end
C
C      subroutine xyl(nlm,np,wp,fp,r2,yl,fl)
CC- Yl-projection of function tabulated on an angular mesh
C      implicit none
C      integer nlm,np,ip,ilm
C      double precision fl(nlm),r2(np),fp(np),yl(np,nlm),wp(np)
C      double precision rl
C      procedure(integer) :: ll
C
C      call dpzero(fl,nlm)
C      do  ip = 1, np
C      do  ilm = 1, nlm
C        rl = dsqrt(r2(ip))**ll(ilm)
C        fl(ilm) = fl(ilm) + fp(ip)*wp(ip)*yl(ip,ilm)/rl
C      enddo
C      enddo
C      end
C      subroutine pryl(mode,strn,mlm,pm,fl,nlm,fac)
CC- Prints out nonzero (l,m) components of fl(ilm)*fac
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 print fac*fl
CCi         :1 print fac*sign(fl)*fl**2
CCi         :1 print fac*sign(fl)*fl**2*(2l+1)*l
CCi   strn
CCi   mlm   :ilm from which this function was derived
CCi   pm    :1   delta m = m(ilm) - m(fl)
CCi         :-1  delta m = m(ilm) - m(fl)
CCi         :0   print l and m not delta l and delta m
CCi         :2   print l and -m not delta l and delta m
CCi   fl    :function as linear combation of Ylm
CCi   nlm   :print coffs to fl up to nlm
CCi   fac
CCo Outputs
CCs Command-line switches
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   06 Apr 18
CC ----------------------------------------------------------------------
C      implicit none
C      character*(*) strn
C      integer mode,mlm,nlm,pm
C      double precision fl(nlm),fac
C      character*120 outs
C      integer ilm,l2,m2,lmax,l,m,lav,dl,dm
C      double precision f
C      real(8),parameter:: tol=1d-6
C      procedure(integer) :: ll
C
C      l = ll(mlm)
C      lav = l*l+l+1
C      m = mlm-lav
C      call awrit2(strn,outs,len(outs),0,l,m)
C
C      lmax = ll(nlm)
C      ilm = 0
C      do  l2 = 0, lmax
CC        ilm0 = l2**2 + l2 + 1   ! the m=0 element for this l
CC        ilmm = ilm0 - 2*l2      ! the m=0 element for this l-1
CC        ilmp = ilm0 + 2*(l2+1)  ! the m=0 element for this l+1
C        do  m2 = -l2, l2
C          ilm = ilm+1
C          if (abs(fl(ilm)) > tol) then
C            if (mode == 0) f = fl(ilm)*fac
C            if (mode == 1) f = fl(ilm)**2*fac*dsign(1d0,fl(ilm))
C            dl = l2-l; if (pm==0 .or. pm==2) dl = l2
C            dm = m2-m; if (pm<0) dm=m2+m; if (pm==0) dm = m2; if (pm==2) dm = -m2
C            call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,f)
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l-1))
CC           call awrit3('%a %,3i%,3i%;12,6D',outs,len(outs),0,dl,dm,fl(ilm)**2*2*(2*l+3))
C
CC           if (l2-l /= -1) stop 'oops'  ! True for grad r^l YL, YL = SH (not polynomials)
CC           if (l2-l /= 1) stop 'oops'  ! True for grad r^-l-1 YL, YL = SH (not polynomials)
C
C          endif
C        enddo
C      enddo
C      call info0(1,0,0,trim(outs))
C
C      end
C      subroutine prmx(strn,s,ns,nr,nc)
CC- writes matrix into out file (for debugging)
C      implicit none
C      integer nr,nc,ns,ifi
C      double precision s(ns,nc,2)
C      character*(14) fmt, fmt0, strn*(*), outs*100
C      integer i,j,fopna,i1mach
C      save fmt
C      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
C      ifi = fopna('out',29,0)
C      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
C#elseC
C      call awrit2('%% rows %i cols %i real',' ',100,ifi,nr,nc)
C#endifC
C      do  i = 1, nr
C        write (ifi,fmt) (s(i,j,1),j=1,nc)
C      enddo
C      write(ifi,*)
CC      do  12  i = 1, nr
CC   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
C      call fclose(ifi)
C
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#elseC
C      outs = ' prm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-100,-i1mach(2))
C#endifC
C      read(*,'(a100)') outs
C
C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in prmx'
C#elseC
C      if (outs == 'q') call rx0('quit in prmx')
C#endifC
C      return
C
C      entry prmx0(fmt0)
C      fmt = fmt0
C      end
C      subroutine zprm(strn,icast,s,ns,nr,nc)
CC- Print complex matrix to file out
C      implicit none
C      integer icast,nr,nc,ns,ifi
C      double precision s(2,ns,nc)
C      character*(20) fmt, outs*80, strn*(*)
C      integer i,j,fopna,i1mach
C      character fmt0*(*)
C      save fmt
C      data fmt /'(9f18.11)'/
C
C      outs = ' '
C      if (icast == 1)  outs = ' real'   ! print real part only
C      if (icast == 11) outs = ' symm'
C      if (icast == 2)  outs = ' complex'
C      if (icast == 12) outs = ' herm'
C      ifi = fopna('out',29,0)
C      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#elseC
C      call awrit2('%% rows %i cols %i'//trim(outs),' ',80,ifi,nr,nc)
C#endifC
C      do  i = 1, nr
C        write (ifi,fmt) (s(1,i,j),j=1,nc)
C      enddo
C      if (mod(icast,10) > 1) then
C      write(ifi,*)
C        do  i = 1, nr
C          write (ifi,fmt) (s(2,i,j),j=1,nc)
C        enddo
C      endif
C      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#elseC
C      outs = ' zprm: wrote '//strn
C      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endifC
C      read(*,'(a80)') outs
C
C#ifdefC NOAWRITE
C      if (outs == 'q') stop 'quit in zprm'
C#elseC
C      if (outs == 'q') call rx0('quit in zprm')
C#endifC
C
C      return
C
C      entry zprm0(fmt0)
C      fmt = fmt0
C      end
C
C#endif
