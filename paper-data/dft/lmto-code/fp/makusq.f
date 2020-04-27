      subroutine makusq(mode,s_site,s_spec,s_lat,s_ham,nbas,nsites,isite,nlmax,nrsmo,
     .  ndham,nphimx,ndimh,napw,igapw,nev,nsp,nspc,nspa,isp,iq,q,evec,ppnl,aus,cPkLq)
C- Accumulate coefficients of (u,s) in all augmentation spheres at one k-pt
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pusq1 bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rsma lmxa lmxl kmxt a nr rmt lmxb pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb pusq1 bstrux
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat qlat vol plat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qlat plat cg indxcg jcg cy qlv dlv
Cio    Passed to:  pusq1 bstrux hxpbl ghibl hklbl gklbl hxpgbl ghigbl
Cio                hklgbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  pwmode
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci         :10s digit
Ci         :0 No rotation to spherical harmonics
Ci         :1 rotate from real to spherical harmonics, m=l:-l order
Ci         :2 rotate from real to spherical harmonics, m=-l:l order
Ci         :100s digit
Ci         :1 Also generate cPkLq
Ci   nbas  :number of basis atoms
Ci   nsites:If zero, coefficients are made all sites.
Ci         :If nonzero, coefficients are made just for a subset
Ci         :of sites (see isite); nsites is the number of sites
Ci   isite :sites at which to calculate coefficients; see nsites
Ci   nlmax :1st dimension of aus (maximum nlma over all sites)
Ci   nrsmo :dimensions cPkLq ... max number of radial functions to expand envelopes
Ci   ndham :dimensions aus : must be at least hamiltonian rank
Ci   nphimx:max number of partial waves of a given l channel at any site
Ci   ndimh :dimensions evec
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igapw:G vectors in PW basis, units of qlat (gvlst2.f)
ci   nev   :number of eigenvectors for which to accumulate aus
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nspa  :dimensioning of aus and cPkLq; see Remarks
Ci   isp   :current spin index for collinear spin polarized case
Ci         :isp should be 1 in the noncollinear case
Ci   iq    :qp index (used only to address element in aus; see Remarks)
Ci   q     :Bloch vector
Ci   evec  :eigenvectors for this q
Ci   ppnl  :NMTO pot pars; see potpus.f
Ci         :ppnl(:,l,nsp,ibas), i.e. for l,isp,ibas, the following:
Ci         :(1)  = inverse potential function (not necessarily created)
Ci         :(2)  = normalization of phi (= 1 here for now)
Ci         :(3)  = rmax * log derivative of phi at rmax
Ci         :(4)  = rmax * log derivative of phidot at rmax
Ci         :(5)  = phi at rmax
Ci         :(6)  = phidot at rmax
Ci         :(7)  = normalization of phidot
Ci         :(8)  = normalization <gz|gz> for local orbital
Ci         :(9)  = overlap <phi|gz>
Ci         :(10) = overlap <phidot|gz>
Ci         :(11) = phz at rmax for local orbitals defined, before
Ci                 (phi,phidot) subtracted
Ci         :(12) = dphz (slope of phz) at rmax for local orbitals defined
Co Outputs
Co   aus   :val,slo of w.f. at MT sphere surface added to aus; see Remarks
Co         :or coefficients to (phi,phidot) projection of w.f,
Co         :depending on mode.
Co         *See note on dimensioning of aus and cPklq in Remarks
Co  cPkLq  :(optional; see mode) coefficients to P_kL expansion of evec
Co         *See note on dimensioning of aus and cPklq in Remarks
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cr Remarks
Cr   Makes coefficients for projection of wave function onto
Cr   augmented functions (u,s) which is valid inside the MT spheres.
Cr   u and s are linear combinations of and phi,phidot defined as:
Cr   u has val=1, slo=0 at rmax, s has val=0, slo=1
Cr
Cr   For example, for EELS matrix elements <nk|r|core> we will need
Cr    |nk> = \sum_L(au_nkL*u_l*Y_L + as_nkL*s_l*Y_L)
Cr
Cr   These are generated from the potential later (see vcdmel)
Cr   makusq returns the au_nkL and as_nkL at one spin and k-pt for
Cr   each of the sites in the unit cell, but ...
Cr   if nsites=nbas is passed then coeffs made for each site, otherwise
Cr   coeffs are made just for the nsites sites listed in isite (suclst)
Cr
Cr   Rotation to spherical harmonics:
Cr   Rlm = Ylmc for m>0, and Ylms for m<0.
Cr   Sequence of real harmonics Rlm is:  -ls,...-1s,...,0,1c,...,lc
Cr   Order of spherical harmonics Ylm is m=-l..l.  Relation:
Cr   See also questaal.org/docs/numerics/spherical_harmonics/
Cr    Yl-m =        (Rlm - i*Rl-m)/sqrt(2)
Cr    Ylm  = (-1)^m*(Rlm + i*Rl-m)/sqrt(2) = (-1)^m*conjg(Yl-m)
Cr
Cr  *Note on dimensioning of aus and cPklq
Cr   aus is dimensioned aus(nlmax,ndham*nspc,3,nspa,nsites,iq)
Cr   third argument is index to kind of partial wave
Cr   cPkLq is dimensioned cPkLq(nlmax,ndham*nspc,nrsmo,nspa,nsites,iq)
Cr   intended as k index, but the caller can choose iq to any be index
Cr   Typically iq=1   if you need aus,cPkLq for only one qp at a time
Cr             iq=ikp if you need aus,cPkLq for all qp
Cr   nspa is the number of spins for which aus,cPkLq data are stored
Cr   Use nspa=1   to store aus,cPkLq in spin channel 1
Cr       nspa=nsp to store aus,cPKLq in spin channel isp
Cr   You MUST use nspa=2 if nspc=2
Cb Bugs
Cb   Should not contain parameter nkap0
Cu Updates
Cu   06 Nov 18 aus dimensioned with variable nphimx
Cu   14 Jun 18 Revised treatment of spherical harmonics
Cu   05 Jun 14 Updated to include dimensioning parameter nspa
Cu   09 Feb 14 cPkLq generation includes coefficients for head envelope fns
Cu   31 Dec 13 New 100s digit mode to generate cPkLq.  Altered argument list
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   08 Feb 13 Internally shortens q vector
Cu   10 Nov 11 Begin migration to f90 structures
Cu   16 Jul 09 Enable rotation to spherical harmonics
Cu   06 Jan 09 Adapt to include APW basis functions
Cu   08 Jul 08 Dimension aus separately from evec
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   23 Dec 04 Extended to the spin-coupled case
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   21 Aug 03 Restrict to a list of sites (see Remarks and suclst)
Cu   12 Feb 02 Extended to local orbitals
Cu   28 Mar 01 (MvS) Added mode to generate coefficients to phi,phidot
Cu                   Some rearrangement of coefficients.
Cu   19 Feb 01 (MvS) shortened argument list
Cu   21 Nov 00 (ATP) Adapted from fp mkrout
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nsp,nspc,nspa,isp,iq,ndham,ndimh,napw,nphimx,nev,
     .  nlmax,nrsmo,nsites,isite(nsites)
      integer,target :: igapw(3,napw)
      integer,parameter :: n0=10,nppn=12
      double precision q(3),ppnl(nppn,n0,nsp,nbas)
      double complex evec(ndimh*nspc,nev),aus(nlmax,ndham*nspc,nphimx,nspa,nsites,iq),
     .  cPkLq(nlmax,ndham*nspc,nrsmo,nspa,nsites,iq)
C ... Dynamically allocated local arrays
      real(8), allocatable :: rofi(:)
      real(8), allocatable :: fh(:,:,:),xh(:,:,:),vh(:),dh(:)
      real(8), allocatable :: fp(:,:,:),xp(:,:,:),vp(:,:),dp(:,:)
      integer,pointer :: igapwl(:,:)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Local parameters
      integer,parameter :: nkap0=4
      integer lh(nkap0)
C     integer ilm,l,m
      integer i,ib,is,kmax,lmxa,lmxh,lmxl,mode0,mode1,
     .  mode2,nkapi,nr,nkaph,kmxax,ispc,ksp,iv,ii
      double precision eh(n0,nkap0),rsmh(n0,nkap0),rsma,a,rmt,tol,qs(3)
      double complex zv(1)
      parameter (tol=1d-8)
      procedure(integer) :: nglob  ! ll,fopna,rdm
      procedure(real(8)) :: dlength
C ... External calls
      external bstrux,dcopy,fradhd,fradpk,gtbsl1,orbl,pusq1,pusq2,pusq4,
     .         radmsh,rlocb1,rxi,shorigv,shorps,tcn,tcx,uspecb

C ... Setup
C     stdo = lgunit(1)
C     ipr  = iprint()
      call tcn ('makusq')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      if (nspa < nspc) call rx('makusq: dimension mismatch')

C ... Shorten q; shift APW G vectors to correspond
      igapwl => igapw
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        allocate(igapwl(3,napw))
        call shorigv(q,qs,s_lat%plat,napw,igapw,igapwl)
      endif

C --- Start loop over atoms ---
      do  i = 1, nsites
        if (isite(1) == 0) then
          ib = i
        else
          ib = isite(i)
        endif
        is = s_site(ib)%spec
        rsma = s_spec(is)%rsma
        lmxa = s_spec(is)%lmxa
        lmxl = s_spec(is)%lmxl
        kmax = s_spec(is)%kmxt
        a = s_spec(is)%a
        nr = s_spec(is)%nr
        rmt = s_spec(is)%rmt
        call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
        lmxh = s_spec(is)%lmxb
        if (lmxa == -1) cycle

C   --- Set up all radial head and tail functions, and their BC's
        ii = (lmxh+1)*nkapi
        allocate(rofi(nr),vh(ii),dh(ii))
        allocate(fh(nr,0:lmxh,nkap0),xh(nr,0:lmxh,nkap0))
        call radmsh(rmt,a,nr,rofi)
        call fradhd(nkapi,eh,rsmh,lh,lmxh,nr,rofi,fh,xh,vh,dh)

C       ii = (lmxa+1)*(kmax+1)
        allocate(fp(nr,0:lmxa,0:kmax),xp(nr,0:lmxa,0:kmax))
        allocate(vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax))
        call fradpk(kmax,rsma,lmxa,nr,rofi,fp,xp,vp,dp)
C   ... Add to the coefficient for the projection onto (u,s) for this site
        nkaph = nglob('nkaph')  ! number of envelope function types joined to (u,s,phiz)
        if (mode2 > 0) then
          kmxax = nrsmo-1-nkaph  ! Maximum k-cutoff for PkL expansion
          call pusq1(mode0+100*mode2,ib,isp,nspc,nspa,s_ham%iprmb,nkaph,
     .      nlmax,kmxax,lmxh,nbas,s_site,s_spec,s_lat,qs,ndham,ndimh,
     .      napw,igapwl,nphimx,nev,evec,vh,dh,vp,dp,ppnl(1,1,1,ib),
     .      aus(1,1,1,1,i,iq),cPkLq(1,1,1,1,i,iq))

C        Debugging ... print augmented w.f for state ii
C         deallocate(xp); allocate(xp(nr,(1+lmxa)**2,2)); xp=0
C         ii = 6
C         if (ii /= 0) then
C         do  ilm = 1, (1+lmxa)**2
C         l = ll(ilm)
C         do  m = 0, kmax
C           xp(1:nr,ilm,1) = xp(1:nr,ilm,1) +
C     .       fp(1:nr,l,m)*cPkLq(ilm,ii,m,1,i,iq)
C           xp(1:nr,ilm,2) = xp(1:nr,ilm,2) +
C     .       fp(1:nr,l,m)*dimag(cPkLq(ilm,ii,m,1,i,iq))
C         enddo
C         do  m = 1, nkapi
C           xp(1:nr,ilm,1) = xp(1:nr,ilm,1) +
C     .       fh(1:nr,l,m)*cPkLq(ilm,ii,-m,1,i,iq)
C           xp(1:nr,ilm,2) = xp(1:nr,ilm,2) +
C     .       fh(1:nr,l,m)*dimag(cPkLq(ilm,ii,-m,1,i,iq))
C         enddo
C         enddo
C
C         print *, sngl(xp(nr,1:5,1))
C         print *, sngl(xp(nr,1:5,2))
C
CC         call prrmsh('Re sm phi',rofi,xp,nr,nr,(1+lmxa)**2)
CC         call prrmsh('Im sm phi',rofi,xp(1,1,2),nr,nr,(1+lmxa)**2)
C
C         j = fopna('phi',-1,0)
C         if (rdm(j,0,0,' ',fac,nr-1,2+lmxa) < 0)
C     .     call rx('file mismatch')
C         rewind j
C         read(j,*)
C         fp(1,:,:) = 0
C         do k = 2, nr
C           read(j,*) fac, (fp(k,l,0), l=0, lmxa)
C         enddo
C         call fclose(j)
C         j = fopna('phid',-1,0)
C         if (rdm(j,0,0,' ',fac,nr-1,2+lmxa) < 0)
C     .     call rx('file mismatch')
C         rewind j
C         read(j,*)
C         fp(1,:,:) = 0
C         do k = 2, nr
C           read(j,*) fac, (fp(k,l,1), l=0, lmxa)
C         enddo
C         call fclose(j)
C         call dpzero(xp,size(xp))
C         do  ilm = 1, (1+lmxa)**2
C         l = ll(ilm)
C         do  m = 0, 1
C           xp(1:nr,ilm,1) = xp(1:nr,ilm,1) +
C     .       fp(1:nr,l,m)*aus(ilm,ii,m+1,1,i,iq)
C           xp(1:nr,ilm,2) = xp(1:nr,ilm,2) +
C     .       fp(1:nr,l,m)*dimag(aus(ilm,ii,m+1,1,i,iq))
C         enddo
C         enddo
CC         call prrmsh('Re true phi',rofi,xp,nr,nr,(1+lmxa)**2)
CC         call prrmsh('Im true phi',rofi,xp(1,1,2),nr,nr,(1+lmxa)**2)
C         print *, sngl(xp(nr,1:5,1))
C         print *, sngl(xp(nr,1:5,2))
C         print *
C         endif

        else
          call pusq1(mode0,ib,isp,nspc,nspa,s_ham%iprmb,nkaph,nlmax,0,
     .      lmxh,nbas,s_site,s_spec,s_lat,qs,ndham,ndimh,napw,igapwl,nphimx,
     .      nev,evec,vh,dh,vp,dp,ppnl(1,1,1,ib),aus(1,1,1,1,i,iq),zv)
        endif

        deallocate(fp,xp,vp,dp)
        deallocate(rofi,fh,xh,vh,dh)

C   ... debugging
C        call yprmi('aus(val) isp,ib,iq = %3:1i',
C     .    [isp,ib,iq],0,3,aus(1,1,1,isp,ib,iq),0,nlmax,(lmxa+1)**2,nev)
C        call yprmi('aus(slo) isp,ib,iq = %3:1i',
C     .    [isp,ib,iq],0,3,aus(1,1,2,isp,ib,iq),0,nlmax,(lmxa+1)**2,nev)
C        call yprmi('aus(loc) isp,ib,iq = %3:1i',
C     .    [isp,ib,iq],0,3,aus(1,1,3,isp,ib,iq),0,nlmax,(lmxa+1)**2,nev)

C   --- Rotate from real to spherical harmonics ---
C       Order assumed: m,...,-m
C       |Psi>_l = \Sum_{m}(A_l,m * u_l + B_l,m * s_l)*R_l,m --->
C       |Psi>_l = \Sum_{m}(C_l,m * u_l + D_l,m * s_l)*Y_l,m
C       R_l,m and Y_l,m are the real and spherical harmonics respectively.
C
C             | (-1)^m/sqrt(2)*A_l,-m + i*(-1)^m/sqrt(2)*A_l,m , m>0
C       C_l,m=|  A_l,m                                         , m=0
C             |  1/sqrt(2)*A_l,-m     -  i*1/sqrt(2)*A_l,m     , m<0
C
C       The same relationships apply to D and B.
        if (mode1 > 0) then

          do  ispc = 1, nspc
          ksp = min(max(ispc,isp),nspa)
          do  iv = 1, nphimx
C            call yprmi('aus(iphi=%i) before rotation ksp,i,iq = %3:1i',
C     .        iv,[ksp,i,iq],3,aus(1,1,iv,ksp,i,iq),0,nlmax,(lmxa+1)**2,nev)
            call s2sph(11+100*mod(mode1,2),lmxa+1,0,aus(1,1,iv,ksp,i,iq),nlmax,nev,nlmax,nev,aus(1,1,iv,ksp,i,iq))
C            call yprmi('aus(iphi=%i) after rotation ksp,i,iq = %3:1i',
C     .        iv,[ksp,i,iq],3,aus(1,1,iv,ksp,i,iq),0,nlmax,(lmxa+1)**2,nev)
          enddo

          if (mode2 > 0) then
          do  iv = 1, nrsmo
C            call yprmi('cPkLq(iphi=%i) before rotation ksp,i,iq = %3:1i',
C     .        iv,[ksp,i,iq],3,cPkLq(1,1,iv,ksp,i,iq),0,nlmax,(lmxa+1)**2,nev)
            call s2sph(11+100*mod(mode1,2),lmxa+1,0,cPkLq(1,1,iv,ksp,i,iq),nlmax,nev,nlmax,nev,cPkLq(1,1,iv,ksp,i,iq))
C            call yprmi('cPkLq(iphi=%i) after rotation ksp,i,iq = %3:1i',
C     .        iv,[ksp,i,iq],3,cPkLq(1,1,iv,ksp,i,iq),0,nlmax,(lmxa+1)**2,nev)
          enddo
          endif

          enddo


        endif

      enddo

      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        deallocate(igapwl)
      endif

      call tcx('makusq')
      end

      subroutine pusq1(mode,ia,isp,nspc,nspa,iprmb,nkaph,nlmax,kmxax,
     .  lmxh,nbas,s_site,s_spec,s_lat,q,ndham,ndimh,napw,igapw,nphimx,nev,evec,
     .  vh,dh,vp,dp,ppnl,aus,cPkLq)
C- Add to the coefficient for the projection onto (u,s) for one site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Co     Stored:    *
Cio    Passed to: bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa lmxb kmxt rsma rmt
Co     Stored:    *
Cio    Passed to: uspecb bstrux
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: jcg indxcg cy
Co     Stored:    *
Cio    Passed to: bstrux
Ci Inputs
Ci   mode  :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci   ia    :augmentation sphere
Ci   isp   :current spin index for collinear case
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nspa  :the number of spins for which aus,cPkLq data are stored
Ci         :Use nspa=1   to store au,as,az,cPkLq in spin channel 1
Ci         :    nspa=nsp to store au,as,az,cPKLq in spin channel isp
Ci         :You MUST use nspa=2 if nspc=2
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmax :dimensions au,as
Ci   lmxh  :basis l-cutoff
Ci   nbas  :size of basis
Ci   q     :bloch vector
Ci   ndham :dimensions au,as,az
Ci   ndimh :dimension of hamiltonian, evec
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igapw :G vectors in PW basis, units of qlat (gvlst2.f)
Ci   nev   :number of eigenvectors to sum over
Ci   evec  :eigenvectors evec(ndimh,nspc,nev)
Ci         :evecs are supplied for one spin only unless nspc=2
Ci   vh    :value of head function in sphere ia
Ci   dh    :slope of head function in sphere ia
Ci   vp    :value of PkL expansion of tail function in sphere ia
Ci   dp    :slope of PkL expansion of tail function in sphere ia
Ci   ppnl  :NMTO pot pars (potpus.f)
Cl Local variables
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl   kspa  :spin index for storing au,as,sz,cPklq.
Cl         :It is the smaller of ksp,nspa.
Co Outputs
Co   aus   :projection of this evec onto (u,s) functions and possibly LO; see potpus.f
Co         :If mode=1, au = projection of this evec onto (phi,phidot) functions
Cr Remarks
Cr   Adapted from augmbl
Cu Updates
Cu   06 Nov 18 aus dimensioned with variable nphimx
Cu   10 Nov 11 Begin migration to f90 structures
Cu   23 Dec 04 Extended to the spin-coupled case
Cu    4 Jun 04 Relax condition nlmax>=nlma
Cu   10 Apr 02 Redimensionsed eh,rsmh to accommodate larger lmax
Cu   12 Feb 02 Extended to local orbitals
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ia,isp,nspc,nspa,nkaph,lmxh,nlmax,kmxax,
     .  nbas,ndham,ndimh,napw,igapw(3,napw),nphimx,nev,n0,nppn
      parameter (n0=10, nppn=12)
      integer iprmb(ndimh)
      double precision ppnl(nppn,n0,2),vp(*),dp(*),vh(*),dh(*),q(3)
      double complex evec(ndimh,nspc,nev),cPkLq(nlmax,ndham*nspc,-nkaph:kmxax,2),
     .  aus(nlmax,ndham*nspc,nphimx,2)
C     .  as(nlmax,ndham*nspc,3,2),
C     .  az(nlmax,ndham*nspc,3,2)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: b(:),cPkL(:,:),cPhh(:,:)
C#ifdefC ALL3C
C      complex(8), allocatable :: b3(:)
C#endif
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer nkap0,nlmxx
      parameter (nkap0=4,nlmxx=121)
      integer lh(nkap0)
      integer isa,lmxa,lmxha,kmax,nlma,ivec,
     .  ilm,k,ll,nkape,ksp,kspa,ispc,nlmto,mode0,mode2
      double precision rsma,rmt,x,
     .  phi,phip,dphi,dlphi,dphip,dlphip,det,rotp(nlmxx,2,2)
      double precision eh(n0,nkap0),rsmh(n0,nkap0),pa(3)

      call tcn ('pusq1')
      mode0 = mod(mode,10)
      mode2 = mod(mode/100,10)
      isa = s_site(ia)%spec
      pa = s_site(ia)%pos
      lmxa = s_spec(isa)%lmxa
      lmxha = s_spec(isa)%lmxb
      kmax = s_spec(isa)%kmxt
      rsma = s_spec(isa)%rsma
      rmt = s_spec(isa)%rmt
      if (lmxa == -1) return
      nlmto = ndimh-napw

C     nlmha = (lmxha+1)**2
      nlma  = (lmxa+1)**2
C     if (nlma > nlmax) call rxi('makusq: need nlmax',nlma)
C     Count no. envelope functions connecting (phi,phidot) at site ia
      call uspecb(0,1,s_spec,isa,isa,lh,rsmh,eh,nkape)

C --- Make strux to expand all orbitals at site ia ---
      allocate(b((kmax+1)*nlma*ndimh))
      call bstrux(2,s_lat,s_site,s_spec,s_lat%cg,s_lat%indxcg,s_lat%jcg,
     .  s_lat%cy,iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,napw,igapw,b,x)
C#ifdefC ALL3C
C      if (mode2 > 0) then
C        allocate(b3((kmax+1)*nlma*ndimh))
C        call bstrux(12,s_lat,s_site,s_spec,s_lat%cg,s_lat%indxcg,
C     .    s_lat%jcg,s_lat%cy,iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,
C     .    napw,igapw,b3,x)
C      endif
C#endif

C     In noncollinear case, isp=1 always => need internal ispc=1..2
C     ksp is the current spin index in both cases:
C     ksp = isp  in the collinear case
C         = ispc in the noncollinear case
C     whereas ispc=1 for independent spins, and spin index when nspc=2
C     kspa is the smaller of ksp and nspa
C     Thus:
C     ispc is the appropriate index for evec
C     ksp  "  "   "           "     for ppnl
C     kspa "  "   "           "     for au,as,az,cPkLq
      do  ispc = 1, nspc
      ksp = max(ispc,isp)
      kspa = min(ksp,nspa)

      if (mode0 == 1) then
        if (nlma > nlmxx) call rxi('makusq:  nlmxx < nlma=',nlma)
        do  ilm = 1, nlma
          k = ll(ilm)+1
          dlphi  = ppnl(3,k,ksp)/rmt
          dlphip = ppnl(4,k,ksp)/rmt
          phi    = ppnl(5,k,ksp)
          phip   = ppnl(6,k,ksp)
          dphi   = phi*dlphi/rmt
          dphip  = dlphip/rmt*phip
          det    = phi*dphip - dphi*phip
          rotp(ilm,1,1) = dphip/det
          rotp(ilm,1,2) = -dphi/det
          rotp(ilm,2,1) = -phip/det
          rotp(ilm,2,2) = phi/det
        enddo
      endif

C --- Loop over eigenstates ---
      allocate(cPkL(0:kmax,nlma),cPhh(nkaph,nlma))
      do  ivec = 1, nev
C       call pusq3(ndimh,nlma,kmax,evec(1,ispc,ivec),b,cPkL)
        call rlocb1(ndimh,nlma,kmax,evec(1,ispc,ivec),1,b,cPkL)
C       call zprm('cPkL',2,cPkL,kmax+1,kmax+1,nlma)
        call pusq2(mode0,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,nlmax,ndham*nspc,nphimx,
     .    iprmb,cPkL,rotp,evec(1,ispc,ivec),vh,dh,vp,dp,aus(1,ivec,1,kspa))

        if (mode2 > 0) then
C#ifdefC ALL3C
C          cpkl = 0 ! Remake cPkL with envelope expansion fully by PkL
C          call rlocb1(ndimh,nlma,kmax,evec(1,ispc,ivec),1,b3,cPkL)
C#endif
          do  k = 0, kmax
            cPkLq(1:nlma,ivec,k,kspa) = cPkL(k,1:nlma)
          enddo
C#ifndef ALL3C
          cPhh = 0 ! Head part of sm. Hankel expansion
          call pusq4(ia,nkaph,nlmto,nlma,iprmb,evec(1,ispc,ivec),cPhh)
          do  k = 1, nkaph
            cPkLq(1:nlma,ivec,-k,kspa) = cPhh(k,1:nlma)
          enddo
C#endif

        endif

      enddo  ! ivec
      deallocate(cPkL,cPhh)

C  --- Debugging --
C      print *, 'kspa,ib,nev=',kspa,ib,nev
C      call zprm('val',2,au(1,1,1,kspa),nlmax,(lmxa+1)**2,nev)
C      call zprm('slo',2,as(1,1,1,kspa),nlmax,(lmxa+1)**2,nev)
C      call zprm('loc',2,az(1,1,3,kspa),nlmax,(lmxa+1)**2,nev)

      enddo ! ispc

      deallocate(b)
C#ifdefC ALL3C
C     if (mode2 > 0) deallocate(b3)
C#endif
      call tcx('pusq1')

      end

      subroutine pusq2(mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,nlmax,ndhamx,nphimx,
     .  iprmb,cPkL,r,evec,vh,dh,vp,dp,aus)
C- Extract projection of eigenstate onto (u,s,z) for sphere at site ia
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci   ia    :augmentation sphere
Ci   nkape :number of envelope function types which are joined to (u,s)
Ci         :Any ktab > nkape is a local orbital
Ci   kmax  :polynomial cutoff in P_kL expansion of envelope tails
Ci   lmxa  :augmentation l-cutoff
Ci   lmxh  :basis l-cutoff
Ci   nlmto :dimension of lmto component of basis
Ci   nlma  :number of L's in augmentation sphere = (lmxa+1)**2
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus : must be at least hamiltonian rank
Ci   nphimx:dmensions aus: global max # of partial waves of a given l
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co   cPkL  :coefficients to P_kL expansion of evec
Ci   r     :2x2 rotation matrices rotating (phi,phidot) to (u,s)
Ci   evec  :eigenvector
Ci   vh    :value of head function in sphere ia
Ci   dh    :slope of head function in sphere ia
Ci   vp    :value of PkL expansion of tail function in sphere ia
Ci   dp    :slope of PkL expansion of tail function in sphere ia
Co Outputs
Co   aus   :projection of this evec onto (u,s) functions and possibly LO; see potpus.f
Co         :If mode=1, au = projection of this evec onto (phi,phidot) functions
Cl Local variables
Cr Remarks
Cu Updates
Cu   06 Nov 18 aus dimensioned with variable nphimx
Cu   12 Feb 02 Extended to local orbitals
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,nlmax,ndhamx,nphimx,iprmb(*)
      double precision vh(0:lmxh,nkape),dh(0:lmxh,nkape)
      double precision vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
      integer nlmxx
      parameter (nlmxx=121)
      double precision r(nlmxx,2,2)
      double complex aus(nlmax,ndhamx,nphimx),evec(nlmto),cPkL(0:kmax,nlma)
C ... Local parameters
      integer n0,nkap0,norb
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      integer io1,l1,ik1,nlm11,nlm12,ilm1,i1,ilma,k
      integer l,ll
      double precision xx
      double complex wk(nlmxx,2)

C     call tcn('pusq2')
      if (nlmto == 0) return
      if (nlma > nlmxx) call rxi('makusq:  nlmxx < nlma=',nlma)

C --- Loop over all orbitals centered at this site ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

C     Contribution from head part
      do  io1 = 1, norb
        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
C       i1 = hamiltonian offset for first orbital in block
        i1 = offl(io1)-nlm11+1
        if (ik1 <= nkape) then
          do  ilm1 = nlm11, nlm12
            l = ll(ilm1)
            aus(ilm1,1,1) = aus(ilm1,1,1) + vh(l,ik1) * evec(ilm1+i1)
            aus(ilm1,1,2) = aus(ilm1,1,2) + dh(l,ik1) * evec(ilm1+i1)
          enddo
        else
          do  ilm1 = nlm11, nlm12
            aus(ilm1,1,3) = aus(ilm1,1,3) + evec(ilm1+i1)
          enddo
        endif
      enddo

C     Contribution from tail part
      do  ilma = 1, nlma
        l = ll(ilma)
        do  k = 0, kmax
          aus(ilma,1,1) = aus(ilma,1,1) + vp(l,k) * cPkL(k,ilma)
          aus(ilma,1,2) = aus(ilma,1,2) + dp(l,k) * cPkL(k,ilma)
        enddo
      enddo

C     Rotate to (phi,phidot)
      if (mode /= 0) then
        call dcopy(2*nlma,aus(1,1,1),1,wk(1,1),1)
        call dcopy(2*nlma,aus(1,1,2),1,wk(1,2),1)
        do  ilma = 1, nlma
          aus(ilma,1,1) = wk(ilma,1)*r(ilma,1,1) + wk(ilma,2)*r(ilma,2,1)
          aus(ilma,1,2) = wk(ilma,1)*r(ilma,1,2) + wk(ilma,2)*r(ilma,2,2)
        enddo
      endif

C     call tcx('pusq2')
      end

C      subroutine pusq3(ndimh,nlma,kmax,evec,b,a)
CC- Add together coeffs to expand wavefct at this site
C see rlocb1
CCu  Adapted from pvrlakl in nfp, reversing index order in b
C      implicit none
CC ... Passed parameters
C      integer ndimh,nlma,kmax
CC      double complex b(ndimh,nlma,0:kmax),a(0:kmax,nlma),evec(ndimh)
C      double complex b(0:kmax,nlma,ndimh),a(0:kmax,nlma),evec(ndimh)
CC ... Local parameters
C      integer k,ilma,i
C      call dpzero(a, 2*(kmax+1)*nlma)
CC     call tcn('pusq3')
C      do k = 0, kmax
C        do ilma = 1 ,nlma
C          do i = 1, ndimh
CC            a(k,ilma) = a(k,ilma) + evec(i)*b(i,ilma,k)
C            a(k,ilma) = a(k,ilma) + evec(i)*b(k,ilma,i)
C          enddo
C        enddo
C      enddo
CC     call tcx('pusq3')
C      end

      subroutine pusq4(ia,nkaph,nlmto,nlma,iprmb,evec,cPhh)
C- Extract projection of eigenstate onto head envelopes for sphere at site ia
C ----------------------------------------------------------------------
Ci Inputs
Ci   ia    :augmentation sphere
Ci   nkaph :number of envelope function types which are joined to (u,s,phiz)
Ci   nlmto :dimension of lmto component of basis
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   evec  :eigenvector
Co Outputs
Co   cPhh  :Extract portion of this evec with head envelope functions
Cr Remarks
Cr
Cu Updates
Cu   09 Feb 14 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ia,nkaph,nlmto,nlma,iprmb(*)
      double complex evec(nlmto),cPhh(nkaph,nlma)
C ... Local parameters
      integer n0,nkap0,norb
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      integer io1,l1,ik1,nlm11,nlm12,ilm1,i1
      double precision xx

C     call tcn('pusq4')
      if (nlmto == 0) return

C --- Loop over all orbitals centered at this site ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

C     Contribution from head part
      do  io1 = 1, norb
        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
C       i1 = hamiltonian offset for first orbital in block
        i1 = offl(io1)-nlm11+1
        if (ik1 <= nkaph) then
          do  ilm1 = nlm11, nlm12
            cPhh(ik1,ilm1) = evec(ilm1+i1)
          enddo
        endif
      enddo

C     call tcx('pusq4')
      end
