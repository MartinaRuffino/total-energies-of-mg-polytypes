      subroutine smves(mode,nbas,s_site,s_spec,s_lat,k1,k2,k3,qmom,gpot0,vval,
     .  hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,f,rhvsm,zvnsm,zsum,vrmt,qbg)
C- Electrostatic potential of the smooth density.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos class
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  vesgcm mshvmt symvvl ugcomp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  rmt lmxl rg lfoca rfoca qc z ctail etail stc lmxb p
Ci                 pz
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  vesgcm corprm mshvmt symvvl ugcomp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc ng vol alat plat nsgrp qlat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv kv cy symgr ag cg indxcg jcg qlv dlv
Cio    Passed to:  vesft vesgcm mshvmt symvvl ugcomp ggugbl gfigbl
Cio                fklbl gklbl hhugbl hhigbl phhigb hklbl hsmbl hgugbl
Ci Inputs
Ci   mode  :0 use input vconst
Ci         :1 generate vconst as - average v(RMT)
Ci   nbas  :size of basis
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   smrho :smooth density on real-space mesh
Ci   qbg   : back ground charge
Cio Inputs/Outputs
Cio  vconst:constant potential to be added to total
Cio        :On input  vconst is set to a default value
Cio        :On output vconst may be set to the average estat
Cio        :          at the MT boundary.
Co Outputs (see also Remarks)
Co   gpot0 :integrals of compensating gaussians g_RL * phi0~
Co         :For accuracy, integral is split into
Co         :g_RL phi0 (vesgcm) + g_RL [phi0~-phi0] (ugcomp)
Co         :vesgcm projects g_RL to the mesh to do the integral
Co         :ugcomp does its integrals analytically with structure constants
Co         :NB: There is a local analog of gpot0 generated in locpt2.
Co   vval  :coffs to YL expansion of es potential at MT boundary
Co   hpot0 :integrals of semicore smooth Hankels * phi0~
Co   sgp0  :sgp0 = sum_RL integral qmom_RL g_RL phi0~
Co         :     = integral [n0~-n0] phi0~
Co   smpot :smooth potential phi0~ (includes compensating gaussians)
Co   smq   :integral of smooth density n0
Co   qsmc  :pseudocore charge
Co   f     :electrostatic contribution to force.
Co   rhvsm :integral n0~ phi0~ + vconst*smq
Co         :where n0~ is sm. density including compensating gaussians
Co   zvnsm :integral (qcorg-z + rhoc) phi0~
Co   vrmt  :electrostatic potential at rmt, with G=0 term in smpot=0
Cl Local variables
Cl   u00   :integral n0 phi[n0] = n0 phi0
Cl   u0g   :integral n0 [phi0~-phi0)
Cl   ugg   :integral [n0~-n0] [phi0~-phi0]
Cl         :ugg is not used
Cs Command-line switches
Cs   --wsmrho : Write smooth density to file, including compensating gaussians
Cr Remarks
Cr  The total density is a sum of three terms,
Cr
Cr    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
Cr
Cr  The first term is the smooth density on a mesh of points; the
Cr  second is the true density and is defined on a radial mesh for each
Cr  sphere; the last is the 1-center expansion of the smooth density on
Cr  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
Cr  identical to the one-center expansion of n0(mesh).  The sum of the
Cr  three terms converges rapidly with l because errors in n_R(r) are
Cr  mostly canceled by errors in n0_R(r).)
Cr
Cr  We add and subtract a set of compensating gaussian orbitals
Cr
Cr    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
Cr
Cr  which render the integral of the local part (the last 3 terms)
Cr  zero in each RL channel.  The g_RL must be localized enough that
Cr  their spillout beyond the MT radius is negligible.
Cr
Cr  We define
Cr
Cr    n0~ = n0 + compensating gaussians sum_RL Q_RL g_RL
Cr
Cr  In the interstitial, the electrostatic potential of n0~ is the true
Cr  estat potential.  The potential of n0 is called phi0 and the
Cr  potential of n0~ is called phi0~.  The total electrostatic energy
Cr  is computed as
Cr
Cr    the electrostatic energy of  n0~ + integral n0*vconst +
Cr    the electrostatic energy of (neutral) local parts
Cr
Cr  vconst may either be passed as an input (mode=0) or it is
Cr  generated here as the average ves(RMT).
Cr
Cr  This routine computes the estat potential and energy from n0~.
Cr  Some variables used in smves and its subroutines:
Cr    Let n0  = smooth density without the compensating sum_RL Q_RL g_RL
Cr        n0~ = n0 + sum_RL Q_RL g_RL
Cr      phi0  = ves[n0]
Cr      phi0~ = ves[n0~]
Cr      g_RL  = gaussian in RL channel with integral Q_RL
Cr      h_R   = l=0 sm hankel in RL channel, (represents core densities)
Cr    qmom_RL = multipole moment in RL channel of (n_R(r) - n0_R(r))
Cr              so that int n_RL(r)-n0_RL(r) = qmom_RL * g_RL(r)
Cr      gpot0 = vector of integrals g_RL * phi0~
Cr            =  integral g_RL * phi0 [vesgcm]
Cr              +integral g_RL * (phi0~-phi0)    [ugcomp]
Cr               The integral is partitioned to minimize mesh errors.
Cr               The first part is done by projecting g_RL to a mesh
Cr               and integrating the product g_RL*phi0 on the mesh
Cr               The second is done analytically by structure constants
Cr      hpot0 = integral h_R * phi0~ (contributions from core)
Cr            = integral h_R * phi0~  [vesgcm]
Cr             +integral h_R * (phi0~-phi0)    [ugcomp]
Cr       u00   :integral n0 phi0      [vesft]
Cr       u0g   :integral n0 (phi0~-phi0)
Cr       sgp0  :integral [n0~-n0] phi0~
Cr   Therefore :u00 + u0g + sgp0 = integral n0~ [phi0~]
Cr       smq   :integral n0
Cr       vconst:constant potential to be added to total.
Cr             :It is computed from average (v(RMT))
Cr       rhvsm :u00 + u0g + sgp0 + vconst*smq
Cr             := integral [n0 phi0 + n0 [phi0~-phi0] + [n0~-n0] phi0~ + vconst*smq
Cr             := integral n0~ phi0~ + vconst*smq
Cr       zvnsm :integral core density * phi0~ + vconst*(sm-H charge)
Cr  Subroutines called by smves:
Cr    vesft    Adds to the smooth potential the electrostatic contribution from n0
Cr             (i.e. without the compensating gaussians).
Cr             This is very simple; it uses nabla^2 -> G^2 in G-space
Cr
Cr    vesgcm   1. makes the first term in gpot0 (see above)
Cr                = integral g_RL * (phi0 = phi[n0])
Cr             2. makes the first term in hpot0 (see above)
Cr             3. adds [[n0~-n0] to the mesh estat potential
Cr
Cr    ugcomp   1. makes the second term in gpot0
Cr             2. makes the second term in hpot0
Cb Bugs
Cb   It is possible to make vval(l=0) for sites with lmxl=-1, which tells
Cb   value of ves at point.  However, vval doesn't have the
Cb   space allocated, and anyway, vrmt is made for those sites by vesgcm.
Cb   So skip for now.
Cu Updates
Cu   16 Jun 16 New ability to add sm external potential
Cu   25 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxl=-1
Cu   19 Sep 02 (WRL) Added background term
Cu   24 Aug 01 Extended to calc vval.  Altered argument list.
Cu   20 Apr 01 Generates vrmt
Cu   21 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp ves_smooth.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,mode
      double precision qsmc,smq,rhvsm,sgp0,vconst,zsum,zvnsm,qbg
      double precision qmom(*),f(3,nbas),
     .  gpot0(*),vval(*),hpot0(nbas),vrmt(nbas)
      double complex smrho(k1,k2,k3,2),smpot(k1,k2,k3,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      complex(8), allocatable :: cv(:,:,:),cgsum(:)
C ... Local parameters
      integer ib,ilm,ipr,is,iv0,lfoc,lgunit,lmxl,n1,n2,
     .  n3,ng,ngabc(3),nglob,nlm,nsp,stdo,j1,j2,j3
      integer procid
      integer, parameter :: master=0
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision R,ceh,cofg,cofh,eint,hsum,pi,qcorg,qcorh,qsc,
     .  rfoc,rmt,s1,s2,sbar,srfpi,sum1,sum2,u00,u0g,ugg,usm,vbar,
     .  vcnsto,vol,xx,y0,z
      character*80 outs
      procedure(logical) :: cmdopt
      procedure(integer) :: iprint,mpipid
      procedure(complex(8)) :: zdotc

C ... Setup
      call tcn('smves')
      procid = mpipid(1)
      ipr   = iprint()
      stdo  = lgunit(1)
      nsp   = nglob('nsp')
      pi    = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0    = 1d0/srfpi
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol

C --- Smooth electrostic potential, and energy integrals ---
C ... FT of smooth density to reciprocal space
C     Electrostatics depend only on total spin density
      if (nsp == 2) then
        call daxpy(k1*k2*k3*2,1d0,smrho(1,1,1,2),1,smrho,1)
      endif
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)

C ... cv = estatic potential of smooth density without gaussians
      allocate(cv(ng,nsp,2))  ! For now, cv can only be one spin.
      call dpzero(cv,2*ng*nsp*2)
      call vesft(s_lat,ng,s_lat%gv,s_lat%kv,k1,k2,k3,smrho,cv,u00)

C ... Debugging check
C      call gvputf(ng,1,s_lat%kv,k1,k2,k3,cv,smpot)
C      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
C      call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
C      call chk_grfmsh(s_lat,k1,k2,k3,smpot,smrho)

C ... Integrals of compensating gaussians with ves; add ves[n0~-n0] to smpot through cv
C     On exit:
C     gpot0 = integrals [ (compensating gaussians g_RL) * phi(n0)]
C     hpot0 = integrals [ (sm hankels) * phi(n0)]  for core
C     ugcomp (below) adds integrals g_RL * [phi(n0~) - phi(n0)] and h_R * [phi(n0~) - phi(n0)]
C     cv    = ves[n0~]
C     vrmt is made but not used; it is overwritten by symvvl
      allocate(cgsum(ng))
      call dpzero(f,3*nbas)
      call vesgcm(1,s_site,s_spec,s_lat,nbas,qmom,ng,
     .  cv,cgsum,f,gpot0,hpot0,qsmc,zsum,vrmt)

      call info5(40,1,0,' Smooth charges: Qmesh = %;6;6d'//
     .  '  Qgauss = %;6,6d  core-nuc = %;6;6d  tot = %;6;6d',
     .  smrho(1,1,1,1)*vol,dble(cgsum(1)*vol),qsmc-zsum,
     .  dble((smrho(1,1,1,1)+cgsum(1))*vol)+qsmc-zsum,0)

      call prfrce(50,nbas,f,'%N After vesgcomp: forces are:')
C     call yprm('gpot0',1,gpot0,0,s_pot%nlml,s_pot%nlml,1)
C     call zprm3('smpot',0,smpot,n1,n2,n3)

C ... Scatter smooth ves into smpot
      call gvputf(ng,1,s_lat%kv,k1,k2,k3,cv,smpot)

C ... Compute ves at MT boundary
      call mshvmt(nbas,s_site,s_spec,s_lat,ng,s_lat%gv,cv,vval)
      call symvvl(nbas,s_site,s_spec,s_lat,vval,vrmt)

C     call zprm3('smpot',0,smpot,n1,n2,n3)

C --- Make vbar = avg v(RMT) and optionally assign to vconst ---
      vbar = 0
      sbar = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        rmt = s_spec(is)%rmt
        vbar = vbar + rmt**2 * vrmt(ib)
        sbar = sbar + rmt**2
      enddo
      vbar = vbar/sbar
      vcnsto = vconst
      if (mode /= 0) vconst = -vbar

      call info5(20,1,0,
     .  ' Average es pot at rmt = %,6;6d  avg sphere pot = %,6;6d   vconst = %,6;6d',vbar,vcnsto,vconst,0,0)

C ... Adjust vbar, vval, gpot0 by vconst
      iv0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        if (s_spec(is)%lmxl < 0) cycle
        nlm = (s_spec(is)%lmxl+1)**2
        vrmt(ib) = vrmt(ib) + vconst
        vval(1+iv0) = vval(1+iv0) + vconst/y0
        gpot0(1+iv0) = gpot0(1+iv0) + vconst/y0
        iv0 = iv0 + nlm
      enddo
      if (ipr >= 40) then
        write (stdo,233)
  233   format(' Average es pot at MT boundaries after vconst shift')
        call arrprt(' Site    ves','%,4i%:-3,6;6d','Id',
     .    nbas,0,4,0,'  | ',xx,vrmt,xx,xx,xx,xx,xx,xx)
      endif

C ... Add background to smrho
      if (qbg /= 0) then
        smrho(1,1,1,1) = smrho(1,1,1,1) + qbg/vol
        R = (3d0/pi/4d0*vol)**(1d0/3d0)
        eint = qbg*2*9d0/10d0/R
        call info(30,0,0,' cell interaction energy from homogeneous'//
     .    ' background (q=%d) is %;6,6d',qbg,eint)
      endif

C ... Scatter smooth rho with compensating background and printout
      if (procid == master .and. cmdopt('--wsmrho',8,0,outs)) then
        call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smrho,cv)
        call daxpy(2*ng,1d0,cgsum,1,cv,1)
        deallocate(cgsum); allocate(cgsum(k1*k2*k3))
        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cv,cgsum)
        call fftz3(cgsum,n1,n2,n3,k1,k2,k3,1,0,1)
        call zprm3('$%N Dumping smrho to file out%N',0,cgsum,k1,k2,k3)
      endif

      deallocate(cgsum,cv)

C ... Back transform of density and potential to real-space mesh
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
      call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)

C ... Integrals n0 and n0 phi0~]
      call mshint(vol,1,n1,n2,n3,k1,k2,k3,smrho,sum1,sum2)
      smq = sum1
      call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smrho,smpot,s1,s2)
      u0g = s1 - u00 !  n0 [phi0~ - phi0]

C ... Add g_RL * (phi0~-phi0) and h_R * (phi0~-phi0) to gpot0 and hpot0
C     and corresponding gradients into the force
C     On exit gpot0 = g_RL * (phi0~) and hpot0 = h_R * (phi0~)
      call ugcomp(nbas,s_site,s_spec,s_lat,qmom,gpot0,hpot0,ugg,f)

C --- Collect energy terms; make zvnuc for smooth problem ---
      zvnsm = 0d0
      sgp0 = 0d0
      iv0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        lmxl = s_spec(is)%lmxl
        if (lmxl < 0) cycle
        nlm = (lmxl+1)**2
C       hsum = integral of charge in sm. Hankel
        hsum = -srfpi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
        hpot0(ib) = hpot0(ib) + vconst*hsum
        zvnsm = zvnsm + (qcorg-z)*y0*gpot0(iv0+1) + cofh*hpot0(ib)
        do  ilm = 1, nlm
          sgp0 = sgp0 + qmom(iv0+ilm)*gpot0(iv0+ilm)
        enddo
        iv0 = iv0+nlm
      enddo
      rhvsm = u00 + u0g + vconst*smq + sgp0
      usm = 0.5d0*(rhvsm+zvnsm)

      call info5(51,1,0,
     .  ' Decomposition of electrostatic integrals:%N'//
     .  ' n0  phi0            %;13,6D%N'//
     .  ' n0 (phi0~ - phi0)   %;13,6D%N'//
     .  ' (n0~-n0) phi0~      %;13,6D%N'//
     .  ' vconst*smq          %;13,6D',
     .  u00,u0g,sgp0,vconst*smq,0)
      if (ipr >= 45) then
        call info5(45,0,0,
     .  '%?#(n<=50)#%N Electrostatic integrals:#%27f--- sum#%N'//
     .  ' psval*phi0~         %;13,6D%N'//
     .  ' psnuc*phi0~         %;13,6D%N'//
     .  ' sum                 %;13,6D',
     .    ipr,rhvsm,zvnsm,2*usm,0)
      endif
      call info2(30,0,0,' smooth rhoves%;14,6D   charge%;13,6D',usm,smq)
      call prfrce(50,nbas,f,'%N After electrostatics: forces are:')

C ... Subtract background
      do j1 = 1, k1
        do j2 = 1, k2
          do j3 = 1, k3
            smrho(j1,j2,j3,1) = smrho(j1,j2,j3,1) - qbg/vol
          enddo
        enddo
      enddo
      smq = smq-qbg

C ... Restore spin 1 density, copy potential to second spin channel
      if (nsp == 2) then
        call daxpy(k1*k2*k3*2,-1d0,smrho(1,1,1,2),1,smrho,1)
        call dcopy(k1*k2*k3*2,smpot,1,smpot(1,1,1,2),1)
      endif

      call tcx('smves')

      end
      subroutine prfrce(ipr,nbas,f,msg)
C- Printout forces
      implicit none
      integer ipr,nbas
      character*(*) msg
      double precision f(3,nbas)
      integer ib
      procedure(integer) :: iprint

      call info0(ipr,0,0,msg)
      do  ib = 1, nbas
        call info2(ipr,0,0,'%,4i%3;12,6D',ib,f(1,ib))
      enddo
      end
