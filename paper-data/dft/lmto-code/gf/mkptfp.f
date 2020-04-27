      subroutine mkptfp(s_site,s_pot,nl,nbas,nsp,lrel,lso,ipc,ipdc,lihdim,pfdim,
     .  indxsh,iopt,pp,sop,socscl,pfz,vRLshf,vshft,z,P)
C- Potential functions parameterized from pot pars, complex energy, DLM case
C-----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb ncomp dlmcl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr dpfr ddpfr
Cio    Passed to:  mkfrpf mksopf
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pprel sop socscl
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lrel  :2 for Dirac equation
Ci   lso   :.true. for spin-orbit coupling (should not occur with lrel=2)
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   iopt  :a compound of digits determining what is generated:
Ci         10000s digit controls whether and how matrix form P is recorded in s_site(:).
Ci             What pot fun is generated depends on modeP (modeP=10s digit iopt)
Ci             As of now output is always recorded in vector form in array P, for non CPA sites.
Ci             This may change in future.
Ci             For some modeP, P is also recorded in matrix form:
Ci           0 s_site(:)%pfr (modeP=0)  s_site(:)%ddpfr (modeP=2)  s_site(:)%dpfr (modeP=3)
Ci           1 P is recorded in s_site(:)%pfra (modeP=0)
Ci         1000s digit
Ci           0 vRLshf is not used
Ci           1 RL-dependent shift of enu,C supplied by vRLshf
Ci           2 Special mode for DLM (not used for now, reserved for later)
Ci         100s digit
Ci           0 make P for alpha representation (bet=alpha)
Ci           1 make P for bare representation (alpha=0)
Ci           2 make P for gamma representation (alpha=gamma)
Ci           3 make P for gamma representation averaged over spins
Ci          10s digit determines whether to make P or some derivative:
Ci             ... the following is returned, if 1s digit iopt is < 4
Ci           0 P <- potential function P.  May be stored in s_site%pfr (orbital order)
Ci                                                   or: in s_site%pfra (see 10000s digit)
Ci           1 P <- P-dot.                 May be stored in s_site%dpfr (orbital order)
Ci           2 P <- -1/2 P-dotdot/P-dot.   May be stored in s_site%ddpfr (orbital order)
Ci           3 P <- sqrt(P-dot), choosing abs(sqrt(delta))
Ci           4 P <- P^alp/P^bet to convert g^alp -> g^bet, where bet = e.g. gamma
Ci           5 P <- P^alp/P^bet*(bet-alp) to convert g^alp -> g^bet
Ci           6 P <- (-1/2 P-dotdot/P-dot)-dot
Ci           7 P <- gamma-alpha (returned as complex numbers)
Ci           8 delta in gamma representation
Ci          1s digit
Ci           0 second-order P
Ci           1 third-order P
Ci           2 first-order P
Ci           3 gam, del taken from pfz(7:8)
Ci           4 Return not P but a potential parameter (depending on 10s digit), e.g. gamma-alpha
Ci   pp    :potential parameters, real, by class (atomsr.f)
Ci   pfz   :potential parameters, complex, hamiltonian (aka downfolding) order (makpfz.f)
Ci   vRLshf:array of orbital- and spin-dependent potential shifts
Ci   vshft :array of site potential shifts
Ci   z     :complex energy
Co Outputs
Co   P     :some potential function, depending on iopt
Co         :P is stored in hamiltionian (aka downfolding) order for non-CPA sites
Co         :P is also stored for all sites, including CPA sites, in orbital order:
Co         : s_site%pfr   (modeP = 0)
Co         : s_site%dpfr  (modeP = 1)
Co         : s_site%ddpfr (modeP = 2)
Cr Remarks
Cr   Second order potential function P:
Cr       Let [] = [(gamma-alpha)(z-C) + Delta]  where C = C^gam  and  Delta = Delta^gam
Cr       Note: [] is related to the overlap matrix, which may be written as
Cr       1/o^alpha = (C-enu) - Delta/(gamma-alpha)
Cr       P^alpha             = (z-C)/[] = (Delta/(z-C) + gamma-alpha)^-1
Cr       P^alpha^-1          = []/(z-C) = Delta/(z-C) + gamma-alpha
Cr       dot-P^alpha         = Delta/[]**2
Cr     -1/2 dotdot-P/dot-P   = (gamma-alpha)/[]
Cr   Third order: z -> z + p(z-enu)^3.  Let e'=1+3p(z-enu)**2.
Cr       dot-P               = e' Delta/[]**2
Cr     -1/2 dotdot-P/dot-P   =-3p(z-enu)/e'+e'(gamma-alpha)/[]
Cr   For normal (not-DLM) sites, P's are stored in downfolding order
Cr       as in the original mkptfp (DLM sites have storage allocated
Cr       for compatibility but are skipped). For DLM sites, the P's
Cr       are stored consecutively for each DLM angle at the end of the
Cr       array, now in orbital order (NOT downfolding).
Cr
Cb Bugs
Cb   This routine should make only sqrt(Pdot) and be consistent with mksopf.
Cb   Thus, modeP=1 should be eliminated in favor of modeP=3.
Cu Updates
Cu  05 Jun 16 10000s iopt=1 -> write P to s_site(:)%pfra
Cu  16 Jan 16 (MvS) gamma2 representation averages gamma over CPA sites
Cu  09 Jul 14 (Kirill) Simplified; potential parameters now also stored in s_site
Cu  22 Feb 13 Bug fix for gamma=2
Cu  15 May 08 (Kirill) First created from mkptfp to handle DLM case
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nbas,nsp,lrel,ipc(nbas),lihdim,indxsh(*),iopt
      integer pfdim,ipdc(*)
      integer npfz,n0
      logical lso
      parameter (npfz=8,n0=10)
      real(8) pp(6,nl,nsp,*),vshft(nbas),sop(*),socscl(*),vRLshf(lihdim,nsp)
      complex(8) z,P(pfdim,nsp),pfz(npfz,lihdim,nsp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      complex(8),allocatable:: pfall(:,:,:),pfi(:,:,:)
C ... Local parameters
      integer lmr,ib,l,m,ic,isp,ipr,lgunit,modeP,irep,irepx,stdo,ivshf,nlib,iso,recs,ll
      real(8) oalpha,enu,Ca,Cg,gam,dela,delg,pa,pgam,alp,bet,
     .  xx,enum,Cam,Cgm,ahold(n0),cpagam(nl,nsp)
      complex(8) zz,phold(n0),zmc,ep,epp,eppp,delz,gamz,wk

C     double precision V0,G0,sG0
      logical order1,order3,orderi,savpfp

C ... For CPA
      integer norb,ncomp,iorb,icomp,i0,ic0
      integer mxorb,nglob
      logical ldsite

C --- Setup and printout  ---
      call getpr(ipr)
C     ipr=55
      stdo = lgunit(1)
      ivshf  = mod(iopt/1000,10)
      irep   = mod(iopt/100,10)
      modeP  = mod(iopt/10,10)
      mxorb  = nglob('mxorb')
      order1 = mod(iopt,10) == 2
      order3 = mod(iopt,10) == 1
      orderi = mod(iopt,10) == 3
      if (mod(iopt,10) == 4) call rx('4s digit iopt outdated')
C     governs where in s_site to store P
      recs = modeP+1; if (modeP == 0 .and. mod(iopt/10000,10) == 1) recs = -1
      savpfp =  .true. ! Always write to P for now
      iso = 0; if (lso) iso = 1; if (lrel == 2) iso = 2

!     if (ipr > 50 .and. lrel /= 2) then
      if (ipr >= 50) then
        irepx = irep ; if (modeP >= 4 .and. modeP /= 6) irepx = -1
        call awrit6(' mkptfp: '//
     .    '%?#n==1#with spin-orbit coupling ##%-1j'//
     .    '%?#n==2#Relativistic (lms) ##'//
     .    '%?#n==0#P-##%-1j'//
     .    '%?#n==1#Pdot-##%-1j'//
     .    '%?#n==2#Pdotdot/Pdot-##%-1j'//
     .    '%?#n==3#sqrt(Pdot)-##%-1j'//
     .    '%?#n==4#P^alpha/P^gamma##%-1j'//
     .    '%?#n==5#(gamma-alpha)*P^alpha/P^gamma##%-1j'//
     .    '%?#n==6#(Pdotdot/Pdot)-dot-##%-1j'//
     .    '%?#n==7#gamma-alpha##'//
     .    '%?#n==0#alpha##%-1j%?#n==1#bare##%-1j%?#n==2#gamma##%-1j%?#n==3#gamma-bar##,'//
     .    '%?#n==0# 2nd-order,##%-1j'//
     .    '%?#n==1# 3rd-order,##%-1j'//
     .    '%?#n==2# 1st-order,##%-1j'//
     .    '%?#n==3# inf-order,##'//
     .    ' z = (%1;6g,%1;6g):',
     .    ' ',100,stdo,iso,modeP,irepx,mod(iopt,10),z,dimag(z))
        if (nsp == 1) then
          write(stdo,'(''  ib   P_l(ib) ...'')')
        else if (.not. lso) then
          write(stdo,'(''  ib comp spin  P_l(ib) ...'')')
        endif
      endif

      if (lrel == 2) then
        call mkfrpf(s_site,nl,nbas,ipc,ipdc,lihdim,pfdim,indxsh,iopt,pp,s_pot%pprel,vshft,z,P)
        return
      elseif (lso) then
C       call mksopfx(s_site,nl,nbas,ipc,ipdc,ldim,indxsh,iopt,pp,s_pot%sop,s_pot%socscl,vshft,z)
        call mksopf(s_site,s_pot,nl,nbas,ipc,ipdc,pfdim,indxsh,iopt,pp,sop,socscl,vshft,z,P)
        return
      endif

      ep = 1
      epp = 0
      eppp = 0
C --- For each site, make P or related function ---
      allocate(pfall(mxorb,nsp,nbas))
      pfall = 0
      do  ib = 1, nbas
        norb = s_site(ib)%norb ; ncomp = s_site(ib)%ncomp
        nlib = sqrt(norb + 0.1)        ! nl for this site
        i0 = mxorb*(ib-1)              ! Offset to this site in indxsh
        allocate(pfi(norb,nsp,ncomp))  ! Temp array to hold P for this site
        ldsite = ncomp > 1          ! True if ib is a CPA site
        ic0 = s_site(ib)%dlmcl         ! Parent class

C       Returns gamma ppar, or if a CPA site, species-averaged gamma
        call cpagamma(ncomp,ic0,s_pot%dlmwt,nl,nsp,pp,cpagam)
C       If spin-average gamma ... do on the fly
C        if (irep == 3) then
C          ahold(1:nl) = (cpagam(1:nl,1) + cpagam(1:nl,nsp))/2
C          cpagam(1:nl,1) = ahold(1:nl)
C          cpagam(1:nl,2) = ahold(1:nl)
C        endif

        do  icomp = 1, ncomp           ! Loop over CPA components
          ic = ic0 + icomp - 1
          do  isp = 1, nsp
            do  iorb = 1, norb
              l = ll(iorb)
              phold(l+1) = 0
C         ... Extract pp(alpha)
              enu = pp(1,l+1,isp,ic) + vshft(ib)
              Ca  = pp(2,l+1,isp,ic) + vshft(ib)
              dela= pp(3,l+1,isp,ic)**2
              pa  = pp(4,l+1,isp,ic)
              gam = pp(5,l+1,isp,ic)
C             gam2= (pp(5,l+1,1,ic) + pp(5,l+1,nsp,ic))/2
              alp = pp(6,l+1,isp,ic)
C 956         format(1X,7F10.6)

C         ... Make pp(gamma)
              xx = 1 + (Ca-enu)*(gam-alp)/dela
              delg = dela*xx**2
              Cg = enu + (Ca-enu)*xx

C       ...   P may be m-dependent if vRLshf is added
              enum = enu
              Cgm  = Cg
              Cam  = Ca
              if (ivshf == 1) then
                enum = enu + vRLshf(indxsh(i0+iorb),isp)
                Cgm  = Cg  + vRLshf(indxsh(i0+iorb),isp)
                Cam  = Ca  + vRLshf(indxsh(i0+iorb),isp)
              elseif (ivshf == 2) then
                call rx(' MKPTFP: not ready for ivshf=2')
              endif

C         ... zz is z or z + p (z - enu)**3
              bet = alp         ! Default screening repsn is alpha
              if (order3) then
                pgam = pa - oalpha(enum,Cam,dela,alp,gam)**2
                zz = z + pgam * (z - enum)**3
                ep = 1 + 3 * pgam * (z - enum)**2
                epp = 6 * pgam * (z - enum)
                eppp = 6 * pgam
              elseif (order1) then
                bet   = cpagam(l+1,isp) ! Shift default screening to gamma
                Cgm   = Cam
                delg = dela
                zz = z
              else
                zz = z
              endif

C        ...  orderi ... (c,gam,del) extracted from pfz
              if (orderi) then
                m = indxsh(i0+iorb)
                zmc   = z - pfz(1,m,isp)
                delz  = pfz(7,m,isp)
                gamz  = pfz(8,m,isp)
                ep = 1
                epp = 0
                eppp = 0
              else
                zmc = zz - Cgm
                delz = delg
                gamz = gam
              endif

C             Select screening representation bet if other than alpha
              if (irep == 1) bet = 0.0d0
              if (irep == 2) bet = cpagam(l+1,isp)
              if (irep == 3) bet = (cpagam(l+1,1) + cpagam(l+1,nsp))/2
              ahold(l+1) = bet

C         ... P
              if (modeP == 0) then
                phold(l+1) = zmc/((gamz-bet)*zmc + delz)
C         ... P-dot
              elseif (modeP == 1) then
                phold(l+1) = ep*delz/((gamz-bet)*zmc + delz)**2
C        ... -1/2 P-dot-dot / P-dot
              elseif (modeP == 2) then
                phold(l+1) = ep*(gamz-bet)/((gamz-bet)*zmc+delz) - epp/ep/2
C         ... (-1/2 P-dot-dot / P-dot)-dot
              elseif (modeP == 6) then
                wk = (gamz-bet)/((gamz-bet)*zmc + delz)
                phold(l+1)= epp*wk - (ep*wk)**2-(eppp*ep-epp**2)/ep**2/2
C         ... sqrt(P-dot)
              elseif (modeP == 3) then
                phold(l+1) = sqrt(ep)*sqrt(delz)/((gamz-bet)*zmc +delz)
C         ... P^alp/P^bet
              elseif (modeP == 4) then
                phold(l+1)= ((gamz-bet)*zmc + delz)/((gamz-alp)*zmc + delz)
C         ... P^alp/P^bet*(bet-alp)
              elseif (modeP == 5) then
                phold(l+1)= ((gamz-bet)*zmc + delz)/((gamz-alp)*zmc + delz)*(bet-alp)
C         ... (bet-alpha)
              elseif (modeP == 7) then
                phold(l+1) = bet-alp
C         ... delta(gamma)
              elseif (modeP == 7) then
                call rx('mkptfp check modeP=7')
                phold(l+1)=delg
              else
                call rxi('mkptfp: bad iopt, ',iopt)
              endif

              pfi(iorb,isp,icomp) = phold(l+1)
            enddo
            m = 0
            if (ldsite) m = icomp
            if (ipr >= 50 .and. nsp == 1) then
              call info5(40,0,0,'%,4i%4p%n:2;9F',ib,2*nlib,phold,4,5)
            elseif (ipr >= 50 .and. nsp == 2) then
              call info5(40,0,0,'%2,4i%,5i%n:2;9F',[ib,m],isp,2*nlib,phold,5)
            endif
            if (ipr > 55) then
              write(stdo,335)  (ahold(l), l=1,nlib)
  335         format(1x,'beta',10(f10.5,12x))
            endif

          enddo
        enddo

C   ... Record potential functions in s_site
        if (recs == 0+1) then
          call zcopy(norb*nsp*ncomp,pfi,1,s_site(ib)%pfr,1)
        elseif (recs == -1) then
          call zcopy(norb*nsp*ncomp,pfi,1,s_site(ib)%pfra,1)
        elseif (recs == 1+1) then
          call zcopy(norb*nsp*ncomp,pfi,1,s_site(ib)%dpfr,1)
        elseif (recs == 2+1) then
          call zcopy(norb*nsp*ncomp,pfi,1,s_site(ib)%ddpfr,1)
        endif
        if (.not. ldsite) pfall(1:norb,:,ib) = pfi(:,:,1)
        deallocate(pfi)
      enddo

C ... For non-CPA sites also record P in downfolding order (for now)
      if (savpfp) then
      do  isp = 1, nsp
        lmr = 0
        do  ib = 1, nbas
          norb = s_site(ib)%norb
          xx = 1
          if (s_site(ib)%ncomp /= 1) xx = 0
          do  iorb = 1, norb
            P(indxsh(lmr+iorb),isp) = xx*pfall(iorb,isp,ib)
          enddo
          lmr = lmr + mxorb
        enddo
      enddo
      endif
      deallocate(pfall)
      end

      subroutine basidx(ib,ncomp,idx,ldsite,ib1,ib2)
C- Returns the formal range of "sites" to use (to sort normal and CPA sites)
C-----------------------------------------------------------------------
      integer ib,ncomp,idx,ib1,ib2
      logical ldsite
      if (ncomp > 1) then
C       ib1..ib2 are the CPA slots for given site to fill P for
        ldsite = .true.
        ib1 = idx + 1
        ib2 = idx + ncomp
      else
        ldsite = .false.
        ib1 = ib
        ib2 = ib
      endif
      end
