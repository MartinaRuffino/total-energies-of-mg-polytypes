      subroutine rhomom(opt,ib1,ib2,s_site,s_spec,s_rhat,qmom,vsum)
C- Multipole moments of valence sphere densities
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
Cr     Elements take same form as s_site%rho{12c}.
Cr     Contents may be linked but need not be.
Cr     Passing s_rhat as separate structure enables contents to differ
Ci Inputs
Ci   opt   :switch specifying what is returned in qmom; see qmom
Ci   ib1,ib2: Make multipole moments for sites in range (ib1:ib2)
Ci   nbas  :size of basis
Co Outputs
Co   qmom  :Case opt=0:
Co         :Multipole moments, stored as a single long vector.
Co         := integral r^l Y_L (rho1-rho2) + core spillout
Co         :Case opt=1:
Co         :qmom(1) = net local valence charge
Co         :Case opt=2:
Co         :qmom(1) = net local charge, including nuclear charge
Co   vsum  :sum over all sites ib of difference vs1_ib - vs2_ib
Co         :vs1_ib = integral in sphere of estat potential[true rho]
Co         :vs2_ib = integral in sphere of estat potential[sm rho]
Cr Remarks
Cu Updates
Cu   30 Aug 13 New argument ib1
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Apr 12 Repackaged radial mesh integration quadrature
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Mar 11 New argument opt
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   11 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp rhomom.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,ib1,ib2
      double precision qmom(*),vsum
C ... For structures
!      include 'structures.h'
      type(str_site):: s_site(*)
      type(str_spec):: s_spec(*)
      type(str_rhat):: s_rhat(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:),rwgt(:),h(:),v(:)
C ... Local parameters
      integer ib,ilm,intopt,ipr,iprint,is,j,j1,lfoc,lmxl,nglob,
     .  nlml,nlmx,nr,nsp,stdo
      parameter(nlmx=81)
      double precision a,ceh,cofg,cofh,fpi,qc,qcore,qcorg,qcorh,
     .  qmloc(nlmx),qnuc,qsc,qval,rfoc,rg,rmt,vs1,vs2,y0,z,qvali,qi

C --- Setup ---
      ipr  = iprint()
      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      fpi  = 16d0*datan(1d0)
      y0   = 1d0/dsqrt(fpi)
      nsp  = nglob('nsp')
      if (ipr >= 40 .and. opt == 0) write(stdo,221)
      if (ipr >= 40 .and. opt /= 0) write(stdo,221) 'sphere Q'

C --- Loop over sites ---
      j1 = 1
      vsum = 0d0
      qval = 0
      qcore = 0
      qnuc = 0
      do  ib = ib1, ib2
        is = s_site(ib)%spec
        lmxl = s_spec(is)%lmxl
        z = s_spec(is)%z
        qc = s_spec(is)%qc
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        rg = s_spec(is)%rg
        if (lmxl == -1) cycle
        allocate(rofi(nr),rwgt(nr),h(nr),v(nr))
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        qc = qcorg+qcorh
        nlml = (lmxl+1)**2
        if (nlml > nlmx) call rxi('rhomom: increase nlmx, need',nlml)
        call radmsh(rmt,a,nr,rofi)
        call radwgt(intopt,rmt,a,nr,rwgt)
        call pvrhom(0,nlml,nr,nsp,rofi,rwgt,
     .    s_rhat(ib)%rho1,s_rhat(ib)%rho2,s_rhat(ib)%rhoc,
     .    cofg,cofh,rg,ceh,rfoc,z,qmloc)
        if (opt == 0) call dcopy(nlml,qmloc,1,qmom(j1),1)
        qvali = qmloc(1)/y0
        qi = qvali + qc-z
        qval = qval + qvali
        qcore = qcore + qc
        qnuc = qnuc + z
        if (ipr >= 40) then
          if (opt == 0)
     .    write(stdo,220) ib,1,qmloc(1),qvali,qc,z
          if (opt /= 0)
     .    write(stdo,220) ib,1,qmloc(1),qvali,qc,z,qi
          do  ilm = 2, nlml
C           j = j1+ilm-1
            j = ilm
            if (dabs(qmloc(j)) > 1d-6) write(stdo,220) ib,ilm,qmloc(j)
          enddo
        endif
  220   format(i13,i6,f12.6,f12.6,2f9.2,f12.6)
  221   format(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x,
     .     'Qc',8x,'Z':6x,a)
        call pvrhm2(0,z,a,nr,nlml,nsp,rofi,rwgt,s_rhat(ib)%rho1,
     .    s_rhat(ib)%rho2,s_rhat(ib)%rhoc,qmloc,cofg,cofh,
     .    rg,ceh,rfoc,h,v,vs1,vs2)

        vsum = vsum+vs1-vs2
        j1 = j1+nlml

   10   continue
        deallocate(rofi,rwgt,h,v)
      enddo

      if (opt == 1) qmom(1) = qval
      if (opt == 2) qmom(1) = qval+qcore-qnuc
      end

      subroutine pvrhom(opt,nlml,nr,nsp,rofi,rwgt,rho1,rho2,rhoc,
     .  cofg,cofh,rg,ceh,rfoc,z,qmom)
C- Multipole moments for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 return multipole moments
Ci         :1 return qmom for rho1 only
Ci         :2 return qmom for rho2 only
Ci         :10s digit
Ci         :0 return qmom for both spins
Ci         :1 return qmom spin resolved
Ci   nlml  :L-cutoff for charge
Ci   nr    :number of radial mesh points
Ci   nsp   :number of spins
Ci   rofi  :radial mesh points
Ci   rwgt  :radial integration weights
Ci   rho1  :local true density*r**2, tabulated on a radial mesh
Ci   rho2  :local smoothed density*r**2, tabulated on a radial mesh
Ci   ... the following are used only if 1s digit opt=0
Ci   rhoc  :core density
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   rg    :smoothing radius for compensating gaussians
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   z     :nuclear charge
Co Outputs
Co   qmom  :multipole moments for one site
Cr Remarks
Cr   Q_L = integral r^l (rho1-rho2) + l=0 contr. from core spillout
Cr   The core spillout term is:
Cr      qcore(rhoc)-z  - sm_qcore-sm_qnuc
Cu Updates
Cu   22 Nov 09 New argument opt
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nlml,nr,nsp
      double precision ceh,cofg,cofh,rfoc,rg,z
      double precision rofi(*),rwgt(*),qmom(nlml,nsp),rhoc(nr,nsp),
     .  rho1(nr,nlml,nsp),rho2(nr,nlml,nsp)
C ... Local parameters
      integer n0,i,ilm,l,m,lmxl,ll,isp,opt0,opt1,ksp
      parameter (n0=10)
      double precision ag,delq,fac,gnu,pi,qcor1,qcor2,qnuc2,r,ddot,
     .  rhochs,rhocsm,rhonsm,srfpi,sum(nsp),sumg,xi(0:n0),qcor1s,h(nr)
C     double precision sumh,samh,y0

      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      if (nsp == 1) opt1 = 0
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      lmxl = ll(nlml)
      do  i = 1, nr
        h(i) = rwgt(i)
      enddo
      ilm = 0
      do  l = 0, lmxl
        do m = -l, l
          ilm = ilm+1
          do  isp = 1, nsp
            if (opt0 == 0) then
              sum(isp) = ddot(nr,rho1(1,ilm,isp),1,h,1)
     .                 - ddot(nr,rho2(1,ilm,isp),1,h,1)
            elseif (opt0 == 1) then
              sum(isp) = ddot(nr,rho1(1,ilm,isp),1,h,1)
            elseif (opt0 == 2) then
              sum(isp) = ddot(nr,rho2(1,ilm,isp),1,h,1)
            endif
          enddo
          qmom(ilm,1) = 0
          if (opt1 == 1) qmom(ilm,2) = 0
          do  isp = 1, nsp
            ksp = 1
            if (opt1 == 1) ksp = isp
            qmom(ilm,ksp) = qmom(ilm,ksp) + sum(isp)
          enddo
        enddo
        do  i = 1, nr
          h(i) = h(i)*rofi(i)
        enddo
      enddo

      if (opt0 /= 0) return

C ... l=0 includes core; some might be spilled out
      ag = 1d0/rg
      fac = 4*pi*(ag*ag/pi)**1.5d0

C ... Renormalize gaussian
      sumg = 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        sumg = sumg + rwgt(i)*gnu
      enddo
      fac = fac/sumg

C     Make:
C     qcor1 = true core charge inside rmax
C     qcor2 = smooth core charge inside rmax * fac
C     qnuc2 = smooth nuclear charge inside rmax * fac
C     delq = (true core q-z) - (sm core q - sm z) * fac
C     sumh = 0d0
      qcor2 = 0d0
      qnuc2 = 0d0
      qcor1 = 0d0
      qcor1s= 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac * r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,1d0/rfoc,xi,1)
        rhonsm = -z*gnu
        rhochs = srfpi*cofh*xi(0)*r*r
        rhocsm = srfpi*cofg*gnu + rhochs
C       sumh = sumh + rwgt(i)*rhochs
        qcor2 = qcor2 + rwgt(i)*rhocsm
        qnuc2 = qnuc2 + rwgt(i)*rhonsm
        qcor1 = qcor1 + rwgt(i)*rhoc(i,1)
        qcor1s= qcor1s + rwgt(i)*rhoc(i,nsp)
      enddo
      if (nsp == 2) qcor1 = qcor1+qcor1s
C     samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
      delq = qcor1-z - qcor2-qnuc2
      if (opt0 == 0) then
        qmom(1,1) = qmom(1,1) + delq/srfpi
      else
        qmom(1,1) = qmom(1,1) + delq/srfpi/2
        qmom(1,2) = qmom(1,2) + delq/srfpi/2
      endif

C      y0 = 1d0/srfpi
C      samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
C      write(stdo,942) samh,sumh,samh-sumh
C  942 format(' integral in smH core:',2f12.6/
C     .   ' Core spill-out charge is',f12.6)
C      write(stdo,821) delq
C  821 format('delq=',f12.6)

      end

      subroutine pvrhm2(opt,z,a,nr,nlml,nsp,rofi,rwgt,rho1,rho2,rhoc,
     .  qmom,cofg,cofh,rg,ceh,rfoc,h,v,vsum1,vsum2)
C- Integral over electrostatic potential, to fix energy origin
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :10s digit
Ci         :0 qmom is not spin resolved
Ci         :1 qmom is spin resolved  NOT IMPLEMENTED
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge
Ci   nsp   :number of spins
Ci   rofi  :radial mesh points
Ci   rwgt  :radial integration weights
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Ci   rhoc  :core density
Ci   qmom  :multipole moments for one site (only l=0 term used)
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci         := y0 * pseudocore charge
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   rg    :smoothing radius for compensating gaussians
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   h     :work array of dimension nr
Ci   v     :work array (holds spherical estat potential - poiss0.f)
Co Outputs
Co   vsum1 :integral in sphere of electrostatic potential for true rho
Co   vsum2 :integral in sphere of electrostatic potential for smooth rho
Co         :including compensating gaussians
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nr,nlml,nsp
      double precision a,ceh,cofg,cofh,rfoc,rg,vsum1,vsum2,z,rofi(*),
     .  rwgt(*),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),h(1),v(nr),qmom(1),
     .  rhoc(nr,nsp)
C ... Local parameters
      integer n0,i
      parameter (n0=10)
      double precision af,ag,b,cg,fac,facs,fpi,gnu,pi,q1,q2,r,srfpi,
     .  vhrho,vsum,y0,xi(0:n0)

      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi
      b = rofi(nr)/(dexp(a*nr-a)-1d0)

C ... True density
      q1 = 0d0
      facs = 1d0/(3-nsp)
      do  i = 1, nr
        h(i) = facs*(srfpi*(rho1(i,1,1)+rho1(i,1,nsp))
     .                    + rhoc(i,1) + rhoc(i,nsp))
        q1 = q1 + rwgt(i)*h(i)
      enddo
      call poiss0(z,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
      vsum1 = 0d0
      do  i = 2, nr
        r = rofi(i)
        vsum1 = vsum1 + 4*pi*rwgt(i)*r*r*(v(i)-2*z/r)
      enddo

C ... Smooth density, including compensating gaussians
      cg = qmom(1) + cofg - y0*z
      ag = 1d0/rg
      af = 1d0/rfoc
      fac = fpi*(ag*ag/pi)**1.5d0
      q2 = 0d0
      facs = 1d0/(3-nsp)
      do  i = 1, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,af,xi,1)
        h(i) = srfpi*(facs*(rho2(i,1,1)+rho2(i,1,nsp))
     .                + cg*gnu + cofh*xi(0)*r*r)
        q2 = q2 + rwgt(i)*h(i)
      enddo
      call poiss0(0d0,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
      vsum2 = 0d0
      do  i = 2, nr
        r = rofi(i)
        vsum2 = vsum2 + 4*pi*rwgt(i)*r*r*v(i)
      enddo

C      write (stdo,400) q1,q2,vsum1,vsum2
C  400 format(' qsph=',2f12.6,'   vsum=',2f12.6)
      end
