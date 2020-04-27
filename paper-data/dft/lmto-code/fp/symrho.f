      subroutine symrho(s_site,s_spec,s_lat,lf,smrho,s_rhat,qbyl,hbyl,f)
C- Symmetrize charge density and related quantities
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  symrat prrhat
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl lmxa nr a rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  symrat prrhat
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc plat qlat nsgrp afmt npgrp ng
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag kv gv ips0 bgv
Cio    Passed to:  symrat symsmr
Cio  s_rhat
Ci     Elts read:  rho1 rho2
Co     Stored:     rho1 rho2 symmetrized on output
Co     Allocated:  *
Cio    Elts passed:rho1 rho2
Cio    Passed to:  symrat psymr1 prrhat
Ci Inputs
Ci   lf    :0 do nothing
Ci         :1s digit  1s bit set : symmetrize smooth density
Ci         :          2s bit set : symmetrize local density
Ci         :          4s bit set : also symmetrize qbyl,hbyl
Ci         :10s digit nonzero    : also symmetrize forces
Ci   smrho :smooth density on uniform mesh
Ci Inputs/Outputs
Cio  smrho :smooth density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cr   Special AFM symmetrization (see also Remarks in mksym.f)
Cr   Performed after normal charge symmetrization if s_lat%nsafm > 0.
Cr   Combines space group op 1 (unit operator) and op s_lat%nsafm.
Cr   A check is made that under this op sites (ib,jb) combined in pairs,
Cr   so that for symmetry purposes every site belongs to a class with 2 elements.
Cr   The two sites are symmetrized with the second atom density spin-flipped.
Cu Updates
Cu   07 Aug 17 Revise AFM symmetrization
Cu   19 Dec 13 Begun nocollinear case (site densities only for now)
Cu   03 Jul 13 Replace f77 pointers with f90 ones
Cu   05 Dec 12 Complete migration to f90 structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   19 Jun 00 Packaged from nfp symrat.f and symsmr.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lf
      double precision f(*),qbyl(*),hbyl(*)
      double complex smrho(*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat) ::  s_lat
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: ipc(:),ipca(:),istab(:,:)
      real(8), allocatable :: g(:,:),ag(:,:)
C ... Local parameters
      integer nsp,nspc,nbas,ngabc(3),n1,n2,n3,k1,k2,k3,nsafm,nrclas
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      procedure(integer) :: nglob,iprint
C ... External calls
      external awrit1,daxpy,dcopy,dpadd,dpzero,dscal,fftz30,gvaddf,
     .         gvgetf,gvsym,info,isanrg,prrhat,psymr0,psymr1,psymr2,
     .         psymrf,psymrq,pxsmr1,pysmr1,radmsh,radwgt,shorho,
     .         sitepack,symprj,symrat,symsmr,tcn,tcx

C ... Setup
      call tcn('symrho')
      call info(30,1,0,' Symmetrize density..',0,0)
      nbas = nglob('nbas'); nsp = nglob('nsp'); nspc = nglob('nspc')
      ngabc = s_lat%nabc
      call fftz30(n1,n2,n3,k1,k2,k3)

C ... Symmetrize local densities
      if (lf > 1) then
        allocate(ipc(nbas))
        call sitepack(s_site,1,nbas,'class',1,ipc,[lf])
        call symrat(s_site,s_spec,nbas,nsp,nspc,lf,
     .    s_lat%plat,s_lat%nsgrp,s_lat%istab,s_lat%symgr,s_lat%ag,
     .    ipc,s_rhat,qbyl,hbyl,f)

C   ... Special AFM symmetry
C       Translation maps jb into ib with spin flip ... average
        nsafm = iabs(s_lat%nsafm)
        if (nsafm /= 0) then
          allocate(ipca(nbas),istab(nbas,2),g(9,2),ag(3,2))
          call suafmsym(s_lat,nbas,ipca,istab,g,ag)
          call symrat(s_site,s_spec,nbas,nsp,nspc,102,
     .      s_lat%plat,2,istab,g,ag,
     .      ipca,s_rhat,qbyl,hbyl,f)
          deallocate(ipca,istab,g,ag)
        endif

!       call pshpr(51)
        if (iprint() > 50) call prrhat(nbas,s_site,s_spec,s_rhat)
        deallocate(ipc)
      endif

C ... Symmetrize interstitial density
      if (mod(lf,2) == 1) then
        call symsmr(s_lat,nsp,nsafm,k1,k2,k3,smrho)
C        if (nsafm /= 0) then
C          call symsmr(s_lat,nsp,nsafm,k1,k2,k3,smrho)
C        endif
      endif

      call tcx('symrho')

      end

      subroutine symrat(s_site,s_spec,nbas,nsp,nspc,lf,
     .  plat,nsgrp,istab,symgr,ag,ipc,s_rhat,qbyl,hbyl,f)
C- Symmetrize the atomic charge densities and the forces.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec class pos
Co     Stored:    *
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxl lmxa nr a rmt
Co     Stored:    *
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: plat qlat nsgrp afmt npgrp istab symgr ag
Co     Stored:    *
Cio    Passed to: *
Cio  s_rhat
Ci     Elts read: rho1 rho2
Co     Stored:    rho1 rho2 symmetrized on output
Cio    Passed to: psymr1
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin channels are coupled, otherwise 1
Ci   lf    :0 do nothing
Ci         :1s digit  2s bit set : symmetrize local density
Ci         :          4s bit set : also symmetrize qbyl,hbyl
Ci         :10s digit nonzero    : also symmetrize forces
Ci         :100s digit nonzero   : special AFM mode
Ci Inputs/Outputs
Cio  s_rhat:vector of offsets containing site density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cu Updates
Cu   05 Dec 12 Complete migration to f90 structures
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,nspc,n0,lf,nsgrp
      parameter (n0=10)
      integer ipc(nbas),istab(nbas,*)
      double precision plat(3,3),symgr(3,3,nsgrp),ag(3,nsgrp)
      double precision f(3,nbas),qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
!      type(str_lat) ::  s_lat
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: ipa(:)
      real(8), allocatable :: pos0(:,:),pos(:),sym(:),rofi(:),rwgt(:)
C ... Local parameters
      integer ib,ib0,ic,intopt,ipr,is,lmxa,lmxl,mxint,nclass,nlml,
     .  nlmx,nr,nrclas,stdo,lf1
      double precision qlat(3,3),xx,a,rmt
      character strn*40
      procedure(integer) :: iprint,ival,nglob
      if (lf < 2) return

      call tcn('symrat')
      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      ipr = iprint()
      lf1 = mod(lf,10)

!     plat = s_lat%plat
      call mkqlat(plat,qlat,xx)
!     qlat = s_lat%qlat
!     ngrp = s_lat%nsgrp
C     afmt = s_lat%afmt
C     npgrp = s_lat%npgrp
!     allocate(ips(nbas),pos0(3,nbas))

C ... Separate sites into symmetry classes
C      allocate(ips(nbas))
C      call symcls(nbas,ngrp,symgr,ag,
C     .  ips,pos0,nclass,ipc)
C      if (ipr >= 40) write(stdo,300) nclass
C  300 format(/' symrat: number of symmetry classes is',i3)

      allocate(pos0(3,nbas))
      call sitepack(s_site,1,nbas,'pos',3,xx,pos0)
      nclass = mxint(nbas,ipc)

C --- Special AFM symmetry ---
C      if (ddot(3,afmt,1,afmt,1) /= 0) then
C      do  ib = 1, nbas
CC       Translation maps jb into ib with spin flip ... average
C        jb = ival(istab,ib+npgrp*nbas)
CC       if (ib >= jb) cycle
C        is = s_site(jb)%spec
C        lmxl0 = s_spec(is)%lmxl
C        ib0 = s_spec(is)%lmxa
C        nr0 = s_spec(is)%nr
C        is = s_site(ib)%spec
C        lmxl = s_spec(is)%lmxl
C        ib0 = s_spec(is)%lmxa
C        nr = s_spec(is)%nr
C        call sanrg(.true.,lmxl,lmxl0,lmxl0,'symrho','lmxl')
C        call sanrg(.true.,nr,nr0,nr0,'symrho','nr')
C        nlml = (lmxl+1)**2
C
C        call psymr2(nr,nsp*nspc,nlml,s_rhat(ib)%rho1,s_rhat(jb)%rho1)
C        call psymr2(nr,nsp*nspc,nlml,s_rhat(ib)%rho2,s_rhat(jb)%rho2)
C
C      enddo
C      endif

C --- Loop over classes ---
      allocate(ipa(nbas),pos(3*nbas))
      do  ic = 1, nclass
        call psymr0(-2,ic,nbas,ipc,pos0,pos,ipa,nrclas)
        if (nrclas > 0) then
          ib0 = ipa(1)
          is = s_site(ib0)%spec
          lmxl = s_spec(is)%lmxl
          lmxa = s_spec(is)%lmxa
          nr = s_spec(is)%nr
          nlml = (lmxl+1)**2
          if (ipr >= 40) write(stdo,800) ic,nrclas,nlml
  800     format(/' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

C   ...   Make the projectors; make to at least to l=1 for forces
          nlmx = max0(nlml,4)
          allocate(sym(nlmx*nlmx*nrclas))
          call symprj(nrclas,nlmx,nsgrp,nbas,istab,symgr,ag,plat,qlat,pos,sym)

C   ...   Apply the projectors to rhoat
          if (lmxl > -1) then
C           Symmetrize local density
            if (mod(lf1/2,2) == 1) then
              call psymr1(nrclas,mod(lf/100,2),ipa,nr,nlml,nsp*nspc,nlmx,sym,s_rhat,1)
              call psymr1(nrclas,mod(lf/100,2),ipa,nr,nlml,nsp*nspc,nlmx,sym,s_rhat,2)
            endif

C           Symmetrize site charges and eigval sum
            if (mod(lf1/4,2) == 1) then ! symmetrize local density
              call psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)
            endif
          endif

C   ...   Symmetrize the forces
          if (mod(lf,100)/10 /= 0) call psymrf(nrclas,ipa,nlmx,sym,f)
          deallocate(sym)
        endif
      enddo
      deallocate(ipa,pos)

      if (ipr >= 80) then
        do  ib = 1, nbas
          is = s_site(ib)%spec
          lmxl = s_spec(is)%lmxl
          lmxa = s_spec(is)%lmxa
          if (lmxa == -1) cycle
          nr = s_spec(is)%nr
          nlml = (lmxl+1)**2
          a = s_spec(is)%a
          nr = s_spec(is)%nr
          rmt = s_spec(is)%rmt
          allocate(rofi(nr),rwgt(nr))
          call radmsh(rmt,a,nr,rofi)
          call radwgt(intopt,rmt,a,nr,rwgt)

          call awrit1(' symm loc output density, ib=%i',strn,80,0,ib)
          call shorho(3,strn,nr,nlml,nsp,rofi,rwgt,s_rhat(ib)%rho1,s_rhat(ib)%rho2)

          deallocate(rofi,rwgt)
        enddo
      endif

      deallocate(pos0)
      call tcx('symrat')

      end

      subroutine psymrf(nrclas,ipa,nlmx,s,f)
C- Symmetrize forces
      implicit none
C ... Passed parameters
      integer nrclas,nlmx,ipa(nrclas)
      double precision s(nlmx,nlmx,nrclas),f(3,*)
C ... Local parameters
      integer ia,ib
      double precision x(3)

      x(1) = 0d0
      x(2) = 0d0
      x(3) = 0d0
      do  ia = 1, nrclas
        ib = ipa(ia)
        x(1)= x(1)+s(4,4,ia)*f(1,ib)+s(4,2,ia)*f(2,ib)+s(4,3,ia)*f(3,ib)
        x(2)= x(2)+s(2,4,ia)*f(1,ib)+s(2,2,ia)*f(2,ib)+s(2,3,ia)*f(3,ib)
        x(3)= x(3)+s(3,4,ia)*f(1,ib)+s(3,2,ia)*f(2,ib)+s(3,3,ia)*f(3,ib)
      enddo
      do  ia = 1, nrclas
        ib = ipa(ia)
        f(1,ib) = (s(4,4,ia)*x(1)+s(2,4,ia)*x(2)+s(3,4,ia)*x(3))*nrclas
        f(2,ib) = (s(4,2,ia)*x(1)+s(2,2,ia)*x(2)+s(3,2,ia)*x(3))*nrclas
        f(3,ib) = (s(4,3,ia)*x(1)+s(2,3,ia)*x(2)+s(3,3,ia)*x(3))*nrclas
      enddo
      end

      subroutine psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)
C- Symmetrize l-decomposed site charges and eval sums
      implicit none
C ... Passed parameters
      integer nrclas,nsp,lmxa,ipa(nrclas),n0
      parameter (n0=10)
      double precision qbyl(n0,nsp,*),hbyl(n0,nsp,*)
C ... Local parameters
      integer stdo,ia,ib,iprint,l,lgunit,isp
      double precision qsum(n0,2),hsum(n0,2),fac

      stdo = lgunit(1)
      call dpzero(qsum,2*n0)
      call dpzero(hsum,2*n0)
      fac = 1d0/nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
        do  l = 0, lmxa
          qsum(l+1,isp) = qsum(l+1,isp) + qbyl(l+1,isp,ib)*fac
          hsum(l+1,isp) = hsum(l+1,isp) + hbyl(l+1,isp,ib)*fac
        enddo
        enddo
      enddo

      if (iprint() >= 40) then
        write(stdo,770) (ipa(ia),ia = 1,nrclas)
  770   format(' symmetrized qbyl,hbyl for class containing ib=',20i3)
        if (nsp == 1) write(stdo,780)
     .    (l,qsum(l+1,1),hsum(l+1,1), l=0,lmxa)
        if (nsp == 2) write(stdo,781)
     .    (l,(qsum(l+1,isp),hsum(l+1,isp),isp=1,nsp), l=0,lmxa)
  780   format(i5,2f12.6)
  781   format(i5,2f12.6,'   spin 2',2f12.6)
      endif

      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
        do  l = 0, lmxa
          qbyl(l+1,isp,ib) = qsum(l+1,isp)
          hbyl(l+1,isp,ib) = hsum(l+1,isp)
        enddo
        enddo
      enddo

      end

      subroutine psymrq1(mode,s_spec,s_site,nbas,nl,nsp,dos)
C- Symmetrize l-decomposed site quantity such as orbital moment or dos
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode   :0 => symmetrize to lmxa as given by s_spec
Ci          :1 => symmetrize to lmxa = nl-1
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Ci   dos   :is symmetrized by class
Cu Updates
Cu   26 Sep 14 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nl,nsp
      double precision dos(nl,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated arrays
      integer, allocatable :: ipc(:),ipa(:)
C ... Local parameters
      integer ia,ib,l,isp,lmxa,nclass,mxint,nrclas,ic,is
      double precision qsum(nl,2),fac,xx(1)

      allocate(ipc(nbas),ipa(nbas))
      call sitepack(s_site,1,nbas,'class',1,ipc,xx)
      nclass = mxint(nbas,ipc)

      do  ic = 1, nclass
        call psymr0(-2,-ic,nbas,ipc,xx,xx,ipa,nrclas)
        if (nrclas <= 0) cycle
        ib = ipa(1)
        if (mode == 0) then
          is = s_site(ib)%spec
          lmxa = s_spec(is)%lmxa
          if (lmxa+1 > nl) call rx('psymrq1: wrongly dimensioned dos')
        else
          lmxa = nl-1
        endif
        if (lmxa < 0) cycle

        call dpzero(qsum,nl*nsp)
        fac = 1d0/nrclas
        do  ia = 1, nrclas
          ib = ipa(ia)
          do  isp = 1, nsp
            do  l = 0, lmxa
              qsum(l+1,isp) = qsum(l+1,isp) + dos(l+1,isp,ib)*fac
            enddo
          enddo
        enddo

        do  ia = 1, nrclas
          ib = ipa(ia)
          do  isp = 1, nsp
            do  l = 0, lmxa
              dos(l+1,isp,ib) = qsum(l+1,isp)
            enddo
          enddo
        enddo

      enddo
      deallocate(ipc,ipa)

      end

      subroutine psymr1(nrclas,lafm,ipa,nr,nlml,nspch,nlmx,sym,s_rhat,icmp)
C- Symmetrize density for one class of atoms
C ----------------------------------------------------------------------
Cio Structures
Cio  s_rhat
Ci     Elts read: rho1 rho2
Co     Stored:    rho1 rho2
Cio    Passed to: *
Ci Inputs
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   lafm  :nonzero => swap spins for atoms in 2nd class
Ci   ipa   :list of sites belonging to class ic (psymr0)
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nspch :Number of spin channels
Ci         :1 in nonmagnetic case
Ci         :2 in spin-polarized case when spin channels are decoupled
Ci         :4 in spin-polarized case when spin channels are coupled
Ci   nlmx  :dimensions sym
Ci   sym   :symmetry projectors for each site within this class (symprj)
Ci   icmp  :1 or 2 to symmetrize rho1 or rho2
Co Outputs
Co   See s_rhat
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Nov 12  f90 structures used for rho
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nrclas,nspch,lafm
      integer ipa(nrclas),nlmx,nr,nlml,icmp
      double precision sym(nlmx,nlmx,nrclas)
C ... For structures
!      include 'structures.h'
      type(str_rhat):: s_rhat(*)
C ... Dynamically allocated local arrays
      real(8), pointer :: rhati(:,:),rho(:,:,:)
C ... Local parameters
      integer stdo,ia,ib,iprint,nn,nglob
      double precision wgt

C ... Accumulate symmetrized true density on first site
      allocate(rho(nr,nlml,nspch))
      call dpzero(rho, nr*nlml*nspch)
      do  ia = 1, nrclas
        ib = ipa(ia)
        if (icmp == 1) rhati => s_rhat(ib)%rho1
        if (icmp == 2) rhati => s_rhat(ib)%rho2
        if (lafm /= 0 .and. ia > 1 .and. nspch > 1) then
          call dswap(nr*nlml,rhati(1,1),1,rhati(1,1+nlml),1)
        endif
        call pxsmr1(1d0,nr,nlml,nspch,sym(1,1,ia),rhati,rho,nn)
        if (lafm /= 0 .and. ia > 1 .and. nspch > 1) then
          call dswap(nr*nlml,rhati(1,1),1,rhati(1,1+nlml),1)
        endif
      enddo

C ... Copy to all sites in class
      wgt = nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        if (icmp == 1) rhati => s_rhat(ib)%rho1
        if (icmp == 2) rhati => s_rhat(ib)%rho2
        call dpzero(rhati, nr*nlml*nspch)
        call pysmr1(wgt,nr,nlml,nspch,sym(1,1,ia),rho,rhati,nn)
        if (lafm /= 0 .and. ia > 1 .and. nspch > 1) then
          call dswap(nr*nlml,rhati(1,1),1,rhati(1,1+nlml),1)
        endif
      enddo

      stdo = nglob('stdo')
      if (iprint() >= 40) write(stdo,100) nn,nlml*nlml
  100 format(' psymr: did',i5,'  of',i5)
      deallocate(rho)

      end

      subroutine symsmr(s_lat,nsp,nsafm,k1,k2,k3,smrho)
C- Symmetrize the smooth charge density
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc ng nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:kv gv ips0 bgv
Cio    Passed to:  *
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   k1..k3:dimensions of smrho for smooth crystal density
Ci   nsafm :if >0, use special 2-op AFM setup (see mksym)
Ci         : and symmetrize spin-down of 2
Cio Inputs/Outputs
Cio  smrho :smooth density on uniform mesh is symmetrized
Cr Remarks
Cu Updates
Cu   07 Aug 17  Adapted for AFM symmetry
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3,nsafm
      double complex smrho(k1,k2,k3,nsp)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      complex(8), allocatable :: cv(:,:)
      complex(8), allocatable :: csym(:,:)
C ... Local parameters
      integer n1,n2,n3,ng,ngrp,ngabc(3),isp
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))

      call tcn('symsmr')
      ngabc = s_lat%nabc
      ng = s_lat%ng
      ngrp = s_lat%nsgrp
      if (nsafm > 0) ngrp = 2
      if (ngrp > 1) then
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3,nsp)
        allocate(cv(ng,2),csym(ng,2))
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1) ! rho(r) -> rho(G)
C       call zprm3('smrho before poke',0,smrho,k1,k2,k3)
        do  isp = 1, nsp
          call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smrho(1,1,1,isp),cv)
          call gvsym(ng,s_lat%gv,s_lat%ips0,s_lat%bgv,cv,csym)
          call dpadd(csym,cv,1,ng*2,-1d0)
          call gvaddf(ng,s_lat%kv,k1,k2,k3,csym,smrho(1,1,1,isp))
        enddo
        if (nsafm > 0 .and. nsp == 2) then
          call gvgetf(ng,2,s_lat%kv,k1,k2,k3,smrho,cv)
          call gvsymafm(ng,s_lat%ips0(1+ng),s_lat%bgv(1+ng),cv)
          call gvputf(ng,2,s_lat%kv,k1,k2,k3,cv,smrho)
        endif
C       call zprm3('smrho after poke',0,smrho,k1,k2,k3)
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)  ! V(G) -> V(r)
        deallocate(cv,csym)

C ... Force density to be real and positive
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3,nsp)
C       do  i23 = 1, k2*k3
C       do  i1  = 1, k1
C         smrho(i1,i23,1) = dble(smrho(i1,i23,1))
C      enddo
C      enddo
      else
        call info(30,1,1,' Smooth density not symmetrized (ngrp=1)',0,0)
      endif

      call tcx('symsmr')

      end

C      subroutine smrpos(smrho,k1,k2,k3,n1,n2,n3,nsp)
CC- Make smrho real and positive
C      implicit none
CC ... Passed parameters
C      integer k1,k2,k3,n1,n2,n3,nsp
C      double complex smrho(k1,k2,k3,nsp)
CC ... Local parameters
C      integer stdo,nglob,i1,i2,i3,nneg,isp
C      double precision rmin,xx
C
C      stdo = nglob('stdo')
C      nneg = 0
C      rmin = 999
C      do  isp = 1, nsp
C        do  i3 = 1, n3
C          do  i2 = 1, n2
C            do  i1 = 1, n1
C              xx = dble(smrho(i1,i2,i3,isp))
C              rmin = min(rmin,xx)
C              if (xx < 0) then
C                nneg = nneg+1
C                xx = 1d-8
C              endif
C              smrho(i1,i2,i3,isp) = xx
C            enddo
C          enddo
C        enddo
C      enddo
C
C      if (nneg > 0) write (stdo,1) nneg,rmin
C    1 format(' rhopos (warning): mesh density negative at',i6,
C     .  ' points.  min=',f13.8)
C
C      end

      subroutine rhoqm(smrho,k1,k2,k3,n1,n2,n3,nsp,vol,qsum)
C- Return charge, magnetic moment of smooth density
C ----------------------------------------------------------------------
Ci Inputs
Ci   smrho :smooth density on uniform mesh
Ci   k1..k3:
Ci   n1..n3:
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   vol   :cell volume
Co Outputs
Co   qsum  :qsum(1) = smrho(+) + smrho(-)
Co         :qsum(2) = smrho(+) - smrho(-) (for nsp=2 only)
Cl Local variables
Cl         :
Cr Remarks
Cr   Input smrho is assumed to be (rho1, rho2)
Cr   If instead smrho=(rho1+rho2,rho1-rho2) => qsum(1,2) = q+amom, q-amom
Cu Updates
Cu   13 Dec 08 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,nsp
      double complex smrho(k1,k2,k3,nsp)
      double precision vol,qsum(2)
C ... Local parameters
      integer i,i1,i2,i3
      double precision sumi,q1,fac

      qsum(1) = 0
      qsum(2) = 0
      fac = vol/(n1*n2*n3)
      q1 = 0
      do  i = 1, nsp
        sumi = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          sumi = sumi + dble(smrho(i1,i2,i3,i))
        enddo
        enddo
        enddo
        if (i == 2) qsum(2) = qsum(2) + q1-sumi
        q1 = sumi
        qsum(1) = qsum(1) + sumi
      enddo
      qsum(1) = fac*qsum(1)
      qsum(2) = fac*qsum(2)
C     write(*,333) qsum
C 333 format(' rhoqm : istl charge, moment = ',2f13.7)

      end


      subroutine psymr2(nr,nsp,nlml,rho1,rho2)
C- Average two local densities, but opposite spins
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml :number of L channels for first density
Cio Inputs/Outputs
Cio  rho1  :first density, dimensioned (nr,nlm1,nsp).
Cio  rho2  :second density, dimensioned (nr,nlm2,nsp)
Cio        :On output
Cio        :  rho1(:,:,1) <- rho1(:,:,1)/2 + rho2(:,:,2)/2
Cio        :  rho1(:,:,2) <- rho1(:,:,2)/2 + rho2(:,:,1)/2
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   10 Apr 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,nlml
C     double precision fac1,fac2
      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp)
C ... Local parameters
      integer isp

      if (nlml <= 0 .or. nsp /= 2) return
      do  isp = 1, nsp
        call dscal(nr*nlml,0.5d0,rho1(1,1,isp),1)
        call daxpy(nr*nlml,0.5d0,rho2(1,1,3-isp),1,rho1(1,1,isp),1)
        call dcopy(nr*nlml,rho1(1,1,isp),1,rho2(1,1,3-isp),1)
      enddo

      end
