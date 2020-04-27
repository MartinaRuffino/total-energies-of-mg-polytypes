      subroutine pnunew(nbas,nsp,s_site,s_spec,pmin,pmax,lfrzw,hab,sab,qbyl,hbyl)
C- Makes new boundary conditions pnu for phi,phidot
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pnu pz v0
Co     Stored:     pnu pz
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa rmt idmod mxcst name z a nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   lfrzw :0, float pnu to band CG, provided IDMOD=0 for that channel
Ci         :1 freeze all pnu for all species.
Ci         :  NB: pnu are also frozen for specific species
Ci         :      that have nonzero 4's bit of species->mxcst.
Ci         :10s digit affects floor for low (free-electron like) pnu
Ci         :0  lm version 6 floor
Ci         :1  lm floor, tailored for LDA
Ci         :2  lm floor, tailored for GW
Ci         :100s digit controls spin averaging of pnz
Ci         :0  pnu,pnz are both spin-dependent
Ci         :1  pnz is spin-averaged
Ci   pmin  :lower bound for fractional part of P
Ci   pmax  :upper bound for fractional part of P
Ci   hab   :<u,s | H | u,s> for each pair uu, us, su, ss; see Remarks
Ci   sab   :<u,s | 1 | u,s>
Ci         :NB: hab and sab are no longer used.
Ci   qbyl  :l-decomposed charge
Ci   hbyl  :l-decomposed eigenvalue sum
Cl Local variables
Cl   lfrzv  if T, freeze valence pnu
Cl   lfrzz  if T, freeze local orbital pnz
Co Outputs
Ci   ssite->pnu :are floated to their band CG
Cr Remarks
Cu Updates
Cu   03 Aug 15 New 100 digits lfrzw for further control of PZ
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 Nov 10 Bug fix: pmin correctly applied
Cu   11 Jan 10 Patch phidx call for deep semicore states
Cu   10 May 09 New 10s digit lfrzw
Cu   28 Jun 06 Handles case idmod=3; New constraints pmin,pmax
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu    9 May 02 Added species-specific freezing of pnu
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   20 Sep 01 Some patches to work in pathological cases
Cu   17 Sep 01 When local orbital present, allow semiore pnz to float
Cu   28 Aug 01 Extended to local orbitals.  For now, freeze pnu
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,lfrzw,n0,nab
      parameter(n0=10,nab=9)
      double precision pmin(n0),pmax(n0),
     .  qbyl(n0,nsp,nbas),sab(nab,n0,nsp,nbas),
     .  hbyl(n0,nsp,nbas),hab(nab,n0,nsp,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: g(:),gp(:),rofi(:)
      real(8), allocatable, target :: v0i(:)
      real(8), pointer     :: p_v0(:,:),v0ic(:)
C ... Local parameters
      logical lpz,lfrzv,lfrzz
      integer idmod(n0),stdo,lwronsk,ipr,ib,is,lmxa,l,ipqn,m,isp,nr,nn,mxcst,fltpnz
      double precision pi,rmt,p1,ebar,a,d0l,pfree,pold,ptry,z,
     .  val(5),slo(5),dl,phi,dphi,pnu(n0,2),pnz(n0,2),pnza(n0)
      double precision ez,umegam,phip,dphip,dlphi,dlphip,cz
      double precision pznew,fi(0:10),gi(0:10),xx,dnz,pflor(n0),plow
      character spid*8
      procedure(integer) :: nglob,lgunit

      call info(30,1,0,' Make new boundary conditions for phi,phidot..',0,0)

      lwronsk = nglob('wronsk')
      stdo = lgunit(1)
      call getpr(ipr)
      pi = 4d0*datan(1d0)
C     Return pflor --- a little higher than free-electron values
      call defpq(20+mod(lfrzw/10,10),0d0,n0-1,1,pflor,0)
      fltpnz = mod(lfrzw/100,10)
      val = 0

C --- For each site, do ---
      do  ib = 1, nbas
        is = s_site(ib)%spec
        pnu = s_site(ib)%pnu
        pnz = s_site(ib)%pz
        pnza(1:n0) = pnz(1:n0,1)/2+pnz(1:n0,nsp)/2
        p_v0 => s_site(ib)%v0
        lmxa = s_spec(is)%lmxa
        rmt = s_spec(is)%rmt
        idmod = s_spec(is)%idmod
        mxcst = s_spec(is)%mxcst
        if (lmxa == -1) goto 10
        spid = s_spec(is)%name
        if (mod(mxcst/4,2) /= 0) call ivset(idmod,1,n0,1)

        if (ipr >= 20) write(stdo,320) ib,is,spid
  320   format(/' site',i5,'   species',i4,':',a)
        if (ipr >= 20) write(stdo,311)
        do  l = 0, lmxa
        do  isp = 1, nsp
          m = l+1
          p1 = 2d10
          pznew = 2d10
C         Initially set lfrzv,lfrzz to external constraints.
          lfrzv = (mod(idmod(m),10) /= 0 .and. mod(idmod(m),10) /= 3)
     .            .or. mod(lfrzw,10) /= 0
          lfrzz = lfrzv .or. fltpnz == 2
          lpz = pnz(m,1) /= 0

C         if (dabs(qbyl(m,isp,ib)) > 1d-8 .and. .not. lpz) then
          if (dabs(qbyl(m,isp,ib)) > 1d-8) then
            ebar = hbyl(m,isp,ib)/qbyl(m,isp,ib)

C       ... Log derivative by direct num integration
            is = s_site(ib)%spec
            z = s_spec(is)%z
            a = s_spec(is)%a
            nr = s_spec(is)%nr
            allocate(g(nr*2),rofi(nr*2),v0i(nr))
            call radmsh(rmt,a,nr,rofi)
            call dpscop(p_v0,v0i,nr,1+nr*(isp-1),1,1d0)

            if (mod(idmod(m),10) == 3) then
              val(1) = rmt
              dl = dtan(pi*(0.5d0-mod(pnu(m,isp),10d0)))
              slo(1) = dl + 1d0
              nn = int(mod(pnu(m,1),10d0))-l-1
              allocate(gp(8*nr))
              call phidx(10*lwronsk,z,l,v0i,0d0,0d0,rofi,nr,2,1d-12,ebar,
     .          val,slo,nn,g,gp,phi,dphi,phip,dphip,xx,xx,xx,xx,xx)
              deallocate(gp)
C         ... cz = estimate for energy of orbital with b.c. connecting
C             to Hankel of energy 0
              dlphi = rmt*dphi/phi
              dlphip = rmt*dphip/phip
              umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
              cz = ebar + umegam
              ebar = cz
C         ... estimate for log der. of wf for constant pot of value C
C             dh_l/dr = l*h_l/r - h_l+1, h=g_l/r**(l+1)
C             when cz -> 0,  dl -> -l-1
C             call bessl(cz*rmt**2,m,fi,gi)
C             dl = (l*gi(l) - gi(l+1))/gi(l)
C             p1 = 0.5d0 - datan(dl)/pi
C             val(1) = rmt
C             slo(1) = dl + 1d0
            endif

            call phidx(10*lwronsk+2,z,l,v0i,0d0,0d0,rofi,nr,0,1d-12,ebar,
     .        val,slo,nn,g,xx,phi,dphi,xx,xx,xx,xx,xx,xx,xx)
C           dphip = (slo(2)-phip)/rmt
            if (nn == int(pnu(m,1))-l-1) then
              dl = rmt*slo(1)/val(1) - 1
              p1 = 0.5d0 - datan(dl)/pi
            elseif (ipr >= 10) then
              call info2(10,0,0,' (warning) failed to find proper '
     .          //'node count for l=%i  ebar=%;4d: pnu not calc',l,
     .          ebar)
            endif

C       ... Estimate new pnz for semicore state
            if (lpz .and. int(mod(pnz(m,1),10d0)) < int(pnu(m,1))) then
              if (fltpnz == 1 .and. nsp == 2) then
                allocate(v0ic(nr))
                call dpscop(p_v0,v0ic,nr,1+nr*(isp-1),1,.5d0)
                call dpsadd(v0ic,p_v0,nr,1,1+nr*(2-isp),.5d0)
                dnz = dtan(pi*(0.5d0-mod(pnza(m),10d0)))
              else
                v0ic => v0i
                dnz = dtan(pi*(0.5d0-mod(pnz(m,isp),10d0)))
              endif

              val(1) = rmt
              slo(1) = dnz + 1d0
              nn = int(mod(pnz(m,1),10d0))-l-1
              allocate(gp(8*nr))
              call phidx(10*lwronsk+100,z,l,v0ic,0d0,0d0,rofi,nr,2,1d-12,ez,
     .          val,slo,nn,g,gp,phi,dphi,phip,dphip,xx,xx,xx,xx,xx)
              deallocate(gp)
C             dphip = (slo(2)-phip)/rmt
              dlphi = rmt*dphi/phi
              dlphip = rmt*dphip/phip
              umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
C         ... cz = estimate for energy of orbital with b.c. connecting
C             to Hankel of energy 0 (=> solution for interstitial
C             is constant potential, value C
              cz = ez + umegam
C         ... estimate for log der. of wf for constant pot of value C
C             Maybe better to recalc. val,slo for this cz like above?
              call bessl(cz*rmt**2,m,fi,gi)
C             dh_l/dr = l*h_l/r - h_l+1, h=g_l/r**(l+1)
C             when cz -> 0,  dl -> -l-1
              dl = (l*gi(l) - gi(l+1))/gi(l)
              pznew = 0.5d0 - datan(dl)/pi
              lfrzv = .true.
              if (.not. associated(v0ic,v0i)) deallocate(v0ic)
            else
              lfrzz = .true.
              pznew = 0.5d0 - datan(dble(l))/pi
              ez = 0
            endif

            deallocate(g,rofi,v0i)

C       ... First root of quad equation is the one we want
C           huu = hab(1,m,isp,ib)
C           hus = (hab(2,m,isp,ib) + hab(3,m,isp,ib))/2d0
C           hss = hab(4,m,isp,ib)
C           suu = sab(1,m,isp,ib)
C           sus = (sab(2,m,isp,ib) + sab(3,m,isp,ib))/2d0
C           sss = sab(4,m,isp,ib)
C           a = hss - ebar*sss
C           b = 2d0*(hus - ebar*sus)
C           c = huu - ebar*suu
C           ddd = b*b - 4*a*c
C           if (ddd >= 0) then
C             x1 = (-b-dsqrt(ddd))/(2*a)
CC            q1 = suu+2*x1*sus + x1*x1*sss
CC            h1 = huu+2*x1*hus + x1*x1*hss
C             p1 = 0.5d0-datan(rmt*x1)/pi
C           endif

          endif

C     ... Free-electron value for pnu
          ipqn = pnu(m,isp)
          d0l = l
          pfree = ipqn + 0.5d0 - datan(d0l)/pi
C     ... Floor to use in place of free-electron value
          plow = ipqn - int(pflor(m)) + pflor(m)
          if (pmin(m) > 0 .and. pmin(m) < 1) then
            plow = max(plow,ipqn+pmin(m))
          endif

C     --- Set the new pnu ---
          pold = pnu(m,isp)
          ipqn = pold
          ptry = pold
          if (dabs(p1) < 1d10) ptry = ipqn+p1
          if (.not. lfrzv) then
            pnu(m,isp) = ptry
C       ... Permit pnu no lower than plow
            if (ptry < plow) pnu(m,isp) = plow
C       ... Permit pnu no higher than pmax
            if (pmax(m) > 0 .and. pmax(m) < 1) then
            if (ptry > ipqn+pmax(m)) pnu(m,isp) = ipqn+pmax(m)
            endif
          endif

          if (ipr >= 20 .and. isp == 1) write(stdo,310)
     .      l,idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
          if (ipr >= 20 .and. isp == 2) write(stdo,410)
     .        idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
  310     format(i2,i6,6f12.6)
  410     format(' spn 2',i2,6f12.6)
  311     format(' l  idmod     ql',9x,'ebar',7x,' pold',8x,
     .       'ptry',8x,'pfree',8x,'pnew')

C     --- Set the new pnz ---
          if (lpz) then
            pold = mod(pnz(m,isp),10d0)
            ipqn = pold
            ptry = pold
            pfree = ipqn + 0.5d0 - datan(d0l)/pi
            if (dabs(pznew) < 1d10) ptry = ipqn+pznew
            if (.not. lfrzz) then
              pnz(m,isp) = ptry + (pnz(m,isp)-mod(pnz(m,isp),10d0))
C         ... Permit pnu no lower than pflor
C             d0l = l
              if (ptry < pfree)
     .          pnz(m,isp) = pfree + (pnz(m,isp)-mod(pnz(m,isp),10d0))
            endif

            if (ipr >= 20 .and. isp == 1) write(stdo,520)
     .        l,idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
            if (ipr >= 20 .and. isp == 2) write(stdo,620)
     .        idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
  520     format(i2,i6,'    sc      ',6f12.6)
  620     format(' spn 2',i2,'    sc      ',6f12.6)

          endif

        enddo                   ! spin
        enddo                   ! l

        s_site(ib)%pnu = pnu
        s_site(ib)%pz  = pnz
   10   continue
      enddo                     ! site

      end
