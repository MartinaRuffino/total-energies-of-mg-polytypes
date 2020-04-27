      subroutine rhogkl(ib1,ib2,nsp,mode,r0,s_site,s_spec,s_rhat,kmax,qkl)
C- G_kL expansion of valence sphere densities
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl z qc a nr rmt rg lfoca rfoca ctail etail stc
Ci                 lmxb p pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  *
Ci Inputs
Ci  ib1,ib2: compute expansion coefficents for sites ib1..ib2
Ci   nsp   :1 make qkl for first spin (possibly the only one)
Ci         :2 make qkl combining spins 1 and 2
Ci   mode  : a compound of digits specifying what is to be included
Ci         : in the expansion coefficients
Ci         : 1s   digit = 1 include local density rho1-rho2
Ci         :              2 include local density rho1
Ci         :              3 include local density rho2
Ci         : 10s  digit = 1 include core density rhoc
Ci                        2 include -1 * core density from sm-hankel
Ci                        3 combination 1+2
Ci         : 100s digit = 1 add -1 * nuclear density Y0 Z delta(r)
Ci         : 1000s digit =0 => expand in G_kl, coffs by integration
Ci         :              1 => expand in P_kl, coffs by integration
Ci                        2 => expand in P_kl, coffs by least squares
Ci                             (NOT IMPLEMENTED)
Ci         : 10000s digit 0 => site density from rho1,rho2,rhoc
Ci                        1 => site density from rho1x,rho2x,rhocx
Ci         : 100000s digit 1 => return spin density
Ci    r0   :use r0 for G_kL smoothing parameter, rather than sspec->rg
Ci   kmax  :make expansion coffs to polynomial cutoff kmax
Co Outputs
Co   qkl  :Expansion coefficients, stored as a single long vector.
Co        := integral pkl Y_L integrand
Co        :where integrand is according to mode
Cr Remarks
Cr    In the spin-polarized case, up- and down- spins are combined.
Cu Updates
Cu   05 Sep 15 New option to write spin density
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   11 Nov 09 New argument r0, also 1000s digit mode
Cu   19 Oct 01 Adapted from rhomom.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ib1,ib2,nsp,mode,kmax
      double precision qkl(0:kmax,*),r0
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: h(:),rofi(:),rwgt(:)
C ... Local parameters
      integer ib,ilm,intopt,ipr,iprint,is,j,j1,k,l,lfoc,lmxl,m,nglob,
     .  nlml,nr,stdo
      double precision z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,ceh,rfoc
      double precision df(0:20)
C     real(8), pointer :: rho1(:,:),rho2(:,:),rhoc(:,:)
C     double precision q,xx,ddot

C --- Setup ---
      ipr  = iprint()
      stdo = nglob('stdo')
      call stdfac(20,df)
      if (ipr >= 40) write(stdo,221)
      intopt = 10*nglob('lrquad')

C --- Loop over sites ---
      j1 = 1
      do  ib = ib1, ib2

C        if (mod(mode/10000,10) == 0) then
C          rho1 => s_site(ib)%rho1
C          rho2 => s_site(ib)%rho2
C          rhoc => s_site(ib)%rhoc
C        else
C          rho1 => s_site(ib)%rho1x
C          rho2 => s_site(ib)%rho2x
C          rhoc => s_site(ib)%rhocx
C        endif

        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        rg = s_spec(is)%rg
        if (lmxl == -1) cycle

        allocate(rofi(nr),rwgt(nr))
        if (r0 /= 0) rg = r0
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        qc = qcorg+qcorh
        nlml = (lmxl+1)**2
        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)
        allocate(h(nr*(kmax+1)*(lmxl+1)))

        call pvrgkl(mode,kmax,nlml,nr,nsp,rofi,rwgt,s_rhat(ib)%rho1,
     .   s_rhat(ib)%rho2,s_rhat(ib)%rhoc,h,cofh,rg,ceh,rfoc,z,qkl(0,j1))

CC       Debugging
C        print *, ddot(nr,rwgt,1,s_rhat(ib)%rho2,1)
C        call prrmsh('rho1-rho2',rofi,
C     .    s_rhat(ib)%rho1-s_rhat(ib)%rho2,nr,nr,nlml*nsp)
C        stop
C        call pkl2ro(1001,ib,rg,kmax,nr,nlml,nsp,rofi,rwgt,
C     .    kmax,nlml,xx,qkl(0,j1),s_rhat(ib)%rho1,xx,q)
C        print*, q
CC        call prrmsh('rho2',rofi,s_rhat(ib)%rho2,nr,nr,nlml*nsp)
C        call prrmsh('sm rho',rofi,s_rhat(ib)%rho1,nr,nr,nlml*nsp)

        deallocate(h)
        if (ipr >= 40) then
          write(stdo,222) ib,0,1,(qkl(k,j1), k=0,kmax)
          ilm = 1
          do  l = 1, lmxl
          do  m = -l, l
            ilm = ilm+1
            j = j1+ilm-1
            if (dabs(qkl(0,j))*df(2*l+1) > 1d-6) write(stdo,220) 0,ilm,
     .        (qkl(k,j)*df(2*l+1),k=0,kmax)

          enddo
          enddo
        endif
  222   format(2x,'ib=',i3,i5,i6,10f12.6)
  220   format(9x,i4,i6,f12.6,10f12.6)
  221   format(/' rhogkl:    k   ilm      qkl (2l+1)!! ...')
C#ifdefC DEBUG
C        if (ib == 1) then
C        print *, 'ib=',ib
C        call prtrkl(mode,kmax,rg,nr,nlml,nsp,rofi,rho1,
C     .    rho2,rhoc,qkl(0,j1))
C        endif
C#endif
        j1 = j1+nlml
        deallocate(rofi,rwgt)
      enddo

      end

      subroutine pvrgkl(mode,kmax,nlml,nr,nsp,rofi,rwgt,rho1,rho2,rhoc,
     .  pkl,cofh,rg,ceh,rfoc,z,qkl)
C- Multipole moments for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : a compound of digits specifying what is to be included
Ci         : in the expansion coefficients
Ci         : 1s   digit = 0 rho1, rho2 not included
Ci         :              1 include local density rho1-rho2
Ci         :              2 include local density rho1
Ci         :              3 include local density rho2
Ci         : 10s  digit = 0 no contribution from core densities
Ci         :              1 include core density rhoc
Ci         :              2 include -1 * core density from sm-hankel
Ci         :              3 combination 1+2
Ci         : 100s digit = 1 add -1 * nuclear density Y0 Z delta(r)
Ci         : 1000s digit =0 => expand in G_kl
Ci         :              1 => expand in P_kl
Ci         :              2 => expand in P_kl, coffs by least squares
Ci         :                   Works so far only for l=0
Ci         : 100000s digit 1 => return spin density
Ci   kmax  :k-cutoff for polynomial expansion of radial part
Ci   nlml  :L-cutoff for charge
Ci   nr    :number of radial mesh points
Ci   nsp   :number of spins
Ci   rofi  :radial mesh points
Ci   rwgt  :radial integration weights
Ci   rho1  :local true density*r**2, tabulated on a radial mesh
Ci   rho2  :local smoothed density*r**2, tabulated on a radial mesh
Ci   rhoc  :core density
Ci   pkl   :help array tabulating polynomials Pkl on radial mesh
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   rg    :smoothing radius for compensating gaussians
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   z     :nuclear charge
Co Outputs
Co   qkl  :expansion coefficients for rho
Cr Remarks
Cr   Q_kL = integral p_kl (rho1-rho2) + l=0 contr. from core spillout
Cr   The core spillout term is:
Cr      qcore(rhoc)-z  - sm_qcore-sm_qnuc
Cr   pvrgkl makes this Q_kL when mode=131; partial contr for other modes
Cr   NB: p0l = a**l and scaling factor for k=0 is 4*pi/(a**l * (2l+1)!!)
Cr       => q0l = 4*pi/(2l+1)!! q_l, where q_l is the multipole moment
Cr   Caution:
Cr    pkl expansion seems only good for small r if rsm is small enough
Cr    If rsm is large, orthogonality breaks down.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kmax,nlml,nr,nsp
      double precision ceh,cofh,rfoc,rg,z
      double precision rofi(*),rwgt(*),qkl(0:kmax,nlml),
     .  rhoc(nr,nsp),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),pkl(nr,0:kmax,0:*)
C ... Local parameters
      integer n0,i,ilm,l,m,lmxl,ll,isp,k,mode3,mode5,k1,k2,kdim,keq
      parameter (n0=10)
      double precision ag,fac,y0,xi(0:n0),fpi,factk,dfact,ddot,
     .  df(0:20),wk(nr),smrch,f1,f2
      double precision wk2(0:kmax,0:nlml),rl,dot3
C     double precision gkl(0:kmax,0:nlml)
      double precision normal(kmax+4,kmax+4),rhs(kmax+4)
C     double precision sumh,samh,y0

C     call prrmsh('rho1',rofi,rho1,nr,nr,nlml*nsp)
C     call prrmsh('rho2',rofi,rho2,nr,nr,nlml*nsp)

      fpi  = 16d0*datan(1d0)
      y0   = 1d0/dsqrt(fpi)
      lmxl = ll(nlml)
      mode3 = mod(mode/1000,10)
      mode5 = mod(mode/100000,10)
      call stdfac(20,df)
      if (mode3 == 0 .or. mode3 == 2) then
C        print *, '!!'; kmax=3
        call vecpkl(rofi,rg,nr,kmax,lmxl,nr,kmax,wk,1,pkl,pkl)
C        call prrmsh('pkl',rofi,pkl,nr,nr,kmax+1)
C        call prrmsh('weights',rofi,rwgt,nr,nr,1)
C        if (mode3 == 2) then
C          do  i = 2, nr
C            rwgt(i) = rwgt(i)/rofi(i)**2
C          enddo
C        endif
      else
C     Non-vectorized form ... should be able to integrate with
C     either pkl -> gkl exp, or with gkl -> pkl exp, but
        do  i = 1, nr
C         call radpkl(rofi(i),rg,kmax,lmxl,kmax,wk2)
          call radgkl(rofi(i),rg,kmax,lmxl,kmax,wk2)
          rl = 1
          do  l = 0, lmxl
            pkl(i,0:kmax,l) = wk2(0:kmax,l)*rl
            rl = rl * rofi(i)
          enddo
        enddo
      endif

C ... rho1-rho2 contribution (or rho1, or rho2, depending on mode)
      call dpzero(qkl,nlml*(kmax+1))
      if (mod(mode,10) > 0) then
      if (mod(mode,10) == 1) then
        f1 = 1
        f2 = -1
      elseif (mod(mode,10) == 2) then
        f1 = 1
        f2 = 0
      elseif (mod(mode,10) == 3) then
        f1 = 0
        f2 = 1
      else
        call rx('rhogkl: bad mode')
      endif
      ilm = 0
      do  l = 0, lmxl
        do  m = -l, l
          ilm = ilm+1
          if (mode3 == 0 .or. mode3 == 1) then
            do  k = 0, kmax
              do  i = 1, nr
C               If rg is small enough, these should all integrate to 1
C               call radgkl(rofi(i),rg,kmax,lmxl,kmax,gkl)
C                if (m == -l) qkl(k,ilm) = qkl(k,ilm) +
C     .            rwgt(i)*rofi(i)**(2+l) * gkl(k,l) * pkl(i,k,l)
                do  isp = 1, nsp
                  qkl(k,ilm) = qkl(k,ilm) + rwgt(i) * pkl(i,k,l) *
     .              (f1*rho1(i,ilm,isp) + f2*rho2(i,ilm,isp))
                  if (mode5 == 1 .and. nsp == 2) then
                    f1 = -f1; f2 = -f2
                  endif
                enddo
              enddo
            enddo
          else
            if (l == 0) then
            call dpzero(normal,(kmax+4)*(kmax+4))
            call dpzero(rhs,kmax+4)
            do  k1 = 0, kmax
              do  k2 = 0, kmax
              do  i = 2, nr
                wk(i) = rofi(i)**2 * rwgt(i)
C               Give weight to rim rho more strongly
                wk(i) = rofi(i)**(2+2) * rwgt(i)
              enddo
              normal(k1+1,k2+1) = dot3(nr,rwgt,pkl(1,k1,l),pkl(1,k2,l))
              normal(k1+1,k2+1) = dot3(nr,wk,pkl(1,k1,l),pkl(1,k2,l))
              enddo
              do  isp = 1, nsp
                wk(1) = 0
                do  i = 2, nr
C                  wk(i) =
C     .              (f1*rho1(i,ilm,isp) + f2*rho2(i,ilm,isp))
                  wk(i) = rofi(i)**2 *
     .              (f1*rho1(i,ilm,isp) + f2*rho2(i,ilm,isp))
                enddo
                rhs(k1+1) = rhs(k1+1) + dot3(nr,rwgt,pkl(1,k1,l),wk)
                if (mode5 == 1 .and. nsp == 2) then
                  f1 = -f1; f2 = -f2
                endif
              enddo
            enddo
C           Constraint: fit charge to integrated charge
            rhs(kmax+2) = 0
            do  i = 2, nr
              wk(i) = rofi(i)**2 * rwgt(i)
            enddo
            do  k1 = 0, kmax
              normal(kmax+2,k1+1) = ddot(nr,wk,1,pkl(1,k1,l),1)
              normal(k1+1,kmax+2) = normal(kmax+2,k1+1)
              enddo
              do  isp = 1, nsp
                do  i = 1, nr
                  wk(i) = (f1*rho1(i,ilm,isp) + f2*rho2(i,ilm,isp))
                enddo
                rhs(kmax+2) = rhs(kmax+2) + ddot(nr,rwgt,1,wk,1)
                if (mode5 == 1 .and. nsp == 2) then
                  f1 = -f1; f2 = -f2
                endif
              enddo
C             Constraint: fit rho = true rho at rmax
              rhs(kmax+3) = (f1*rho1(nr,ilm,isp) + f2*rho2(nr,ilm,isp))
              do  k1 = 0, kmax
                normal(kmax+3,k1+1) = pkl(nr,k1,l)*rofi(nr)**2
                normal(k1+1,kmax+3) = pkl(nr,k1,l)*rofi(nr)**2
              enddo
              normal(kmax+3,kmax+3) = 0
C           Constraint: fit slope = slope of true rho at rmax
C            wk = 0
C            do  i = nr-5, nr
C              wk(i) =
C     .         (f1*rho1(i,ilm,isp)+f2*rho2(i,ilm,isp))/rofi(nr)**2
C            enddo
C            call poldif(rofi(nr-5),wk(nr-5),6,4,0,fac,wk)
C            rhs(kmax+4) = wk(6)
C            call vecpkl(rofi(nr),rg,1,kmax,lmxl,1,kmax,wk,11,gkl,wk2)
C            do  k1 = 0, kmax
C              normal(kmax+4,k1+1) = wk2(k1,l)
C              normal(k1+1,kmax+4) = wk2(k1,l)
C            enddo

C           Room for 3 constraints
            kdim = kmax+4
C           Impose 1 constraint
            keq  = kmax+2

            call prmx('normal',normal,kdim,keq,keq)
            call prmx('rhs',rhs,kdim,keq,1)
            call dinv('s',keq,kdim,normal)
C           call prmx('normal^-1',normal,kmax+1,kmax+1,kmax+1)
C           qkl(:,ilm) = matmul(normal,rhs)
            call dgemm('N','N',keq,1,keq,1d0,normal,kdim,
     .        rhs,keq,0d0,qkl(0,ilm),keq)
            call prmx('qkl',qkl(0,ilm),keq,keq,1)
            else
              qkl(:,ilm) = 0
            endif
          endif
        enddo
      enddo
      endif

      if (mode3 == 2) return

C ... Core part (spec'd by 10s digit mode)
      if (mod(mode/10,10) > 0) then
        do  k = 0, kmax
          f1 = 1
          do  isp = 1, nsp
C         Case 1 or 3: add core density
          if (mod(mod(mode/10,10),2) /= 0) then
            do  i = 1, nr
              qkl(k,1) = qkl(k,1) + f1*y0*rwgt(i)*rhoc(i,isp)*pkl(i,k,0)
            enddo
          endif
C         Case 2 or 3: subtract core density from sm. Hankel
          if (mod(mode/10,10) >= 2) then
            do  i = 1, nr
              call hansmr(rofi(i),ceh,1/rfoc,xi,1)
              smrch = cofh*xi(0)*rofi(i)**2
              qkl(k,1) = qkl(k,1) - f1*rwgt(i)*smrch*pkl(i,k,0)
            enddo
          endif
          if (mode5 == 1 .and. nsp == 2) then
            f1 = -f1; f2 = -f2
          endif
          enddo
        enddo
      endif

C ... Nuclear part (spec'd by 100s digit mode)
      if (mod(mode/100,10) == 1 .and. mode5 /= 0) then
        do  k = 0, kmax
          qkl(k,1) = qkl(k,1) - y0*z*pkl(1,k,0)
        enddo
      endif

C ... Scale to get coefficients of the G_kL; see radpkl
      ag = 1/rg
      ilm = 0
      dfact = 1
      do  l = 0, lmxl
        do  m = -l, l
          ilm = ilm+1
          factk = 1d0
          do  k = 0, kmax
            fac = fpi / ((4*ag*ag)**k * ag**l * factk * dfact)
            qkl(k,ilm) = qkl(k,ilm) * fac
            factk = factk*(k+1)
          enddo
        enddo
        dfact = dfact*(2*l+3)
      enddo

      end
C#ifdefC DEBUG
C      subroutine prtrkl(opt,kmax,rg,nr,nlml,nsp,rofi,rho1,rho2,rhoc,qkl)
CC- Printout Gkl or expansion of rho, for debugging
CC  Example for integrating tabulated moment:
CC  mc out.te -e4 x1 x2 'x5*x1' 'x7*x1*x1' -int 0 2.113465
C      implicit none
C      integer opt,kmax,nr,nlml,nsp
C      double precision rg,rofi(nr),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),
C     .  qkl(0:kmax,nlml),rhoc(nr,nsp)
C      integer lmxx,n0
C      parameter(lmxx=6)
C      double precision g(0:kmax,0:lmxx)
C      double precision ,allocatable:: p(:,:,:), rhop(:,:,:)
C      integer lmxl,ll,isp,ilm,ir,k,l,m
C      double precision wk(nr)
C
CC     call prrmsh('rhoc',rofi,rhoc,nr,nr,1)
C      if (mod(opt,10) == 1) then
C        call daxpy(nr*nlml*nsp,-1d0,rho2,1,rho1,1)
C      endif
C      call prrmsh('given rho1-rho2',rofi,rho1,nr,nr,nlml)
C      lmxl = ll(nlml)
C
C      allocate (p(nr,0:kmax,0:lmxl),rhop(nr,nlml,nsp))
C      rhop = 0
C      isp = 1
C      do  ir = 1, nr
C        if (mod(opt/1000,10) == 0) then
C          call radgkl(rofi(ir),rg,kmax,lmxl,kmax,g)
C        else
C          call radpkl(rofi(ir),rg,kmax,lmxl,kmax,g)
C        endif
C        ilm = 0
C        do  l = 0, lmxl
C        do  m = -l, l
C          ilm = ilm+1
C          do  k = 0, kmax
C            rhop(ir,ilm,1) = rhop(ir,ilm,1) +
C     .        qkl(k,ilm) * g(k,l) * rofi(ir)**(2+l)
C          enddo
C        enddo
C        enddo
C      enddo
C      call prrmsh('fit rho',rofi,rhop,nr,nr,nlml)
C      deallocate (p,rhop)
C      if (mod(opt,10) == 1) then
C        call daxpy(nr*nlml*nsp,1d0,rho2,1,rho1,1)
C      endif
C      end
C#endif
