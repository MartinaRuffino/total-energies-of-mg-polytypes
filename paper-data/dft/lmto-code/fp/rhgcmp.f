      subroutine rhgcmp(mode,ib1,ib2,s_site,s_spec,s_lat,s_rhat,kmax,ng,cg)
C- Adds density of compensating gaussians to FT list
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  rhogkl
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl rg lfoca rfoca qc z ctail etail stc lmxb p pz
Ci                 rmt a nr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  corprm rhogkl
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nabc ng vol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:gv
Cio    Passed to:  *
Cio  s_rhat
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rho1 rho2 rhoc
Cio    Passed to:  rhogkl
Ci Inputs
Ci   mode  : a compound of digits specifying what is to be included
Ci         : in the expansion coefficients
Ci         : 1s   digit = 1 add local density rho1-rho2
Ci         :              2 add local density rho1
Ci         :              3 add local density rho2
Ci         : 10s  digit = 1 add core density rhoc
Ci         :              2 add -1 * core density from sm-hankel
Ci         :                in the local density, restoring it
Ci         :                by adding the sm-hankel to the FT mesh
Ci         :              3 combination 1+2
Ci         : 100s digit = 1 add -1 * nuclear density Z delta(r)
Ci         :                In this mode, Z is smoothed into the G_kL
Ci         :              2 add -1 * nuclear density Z delta(r)
Ci         :                In this mode, Z is incporporated directly
Ci         :                in a PW expansion (Z is not smoothed).
Ci         :100000s digit 1 return spin density
Ci         :
Ci         :Examples:
Ci         :mode=130 include the core, the core tail and nuclear charges
Ci         :         This should make the system charge-neutral.
Ci         :mode=131 Like mode=130, but exclude nuclear charge.
Ci         :         The system should have net charge sum_z
Ci         :mode=2   Exclude all core charges, i.e. gaussian (qcorg-z)
Ci         :  and qcorh from the foca Hankel density.
Ci         :  The system should have the valence charge.
Ci         :3 Like 0, but include nuclear charge -Z delta(r)
Ci         :  directly in a PW expansion (Z is not smoothed).
Ci   ng    :number of G-vectors
Co Outputs
Co   cg    :FT of local densities is added to cg, depending on mode.
Cr Remarks
Cr   The local charges inside each augmentation sphere
Cr   (including -1 * the core tail) are smoothed by expanding
Cr   in a  G_kL expansion for k=0..kmax.  The latter is
Cr   subsequently converted into a PW expansion.
Cu Updates
Cu   05 Sep 15 New option to write spin density
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   23 Oct 01 rhgcmp now expands local densities in
Cu             GkL for k=0..kmax, l=1..nlml for each site
Cu             Recovers old rhgcmp for kmax=0.  New argument list.
Cu   09 Feb 01 Added mode
Cu   30 May 00 Adapted from nfp rho_gcomp.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ib1,ib2,ng,kmax
      double complex cg(ng)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat) ::  s_lat
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iv(:)
      real(8), allocatable :: yl(:),g2(:),g(:),cs(:),sn(:),qkl(:)
C ... Local parameters
      integer ib,is,iv0,lmxl,ltop,n1,n2,n3,ng1,nglob,nlm,
     .  nlmtop,nspec,ngabc(3),lfoc,modgkl,nsp
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .  vol,z,q0(3),df(0:20),plat(3,3),qlat(3,3),tau(3)
C ... External calls
      external corprm,poppr,pshpr,rhgcm2,rhgcm3,rhogkl,stdfac,suphas,
     .         suphs0,suylg,tcn,tcx
      data q0 /0d0,0d0,0d0/

      call tcn('rhgcmp')
      call stdfac(20,df)
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      ngabc = s_lat%nabc
      ng1 = s_lat%ng
      vol = s_lat%vol
      nspec = nglob('nspec')
      nsp   = nglob('nsp')
      modgkl = mode
      if (mode >= 200) modgkl = mod(mode,100)
      if (mode >= 100000) modgkl = modgkl + 100000

C      if (mode == 0) then
C        modgkl = 131
C      elseif (mode == 1 .or. mode == 3) then
C        modgkl = 31
C      elseif (mode == 2) then
C        modgkl = 1
C      endif
C      call sanrg(.true.,mode,0,3,'rhgcmp:','mode')

C --- Set up help arrays ---
      ltop = 0
      do  is = 1, nspec
        lmxl = s_spec(is)%lmxl
        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2

      allocate(yl(ng*nlmtop),g2(ng))
      allocate(g(ng*3))
      call suylg(ltop,alat,ng,s_lat%gv,g,g2,yl)
      deallocate(g)

      allocate(iv(ng*3))
      call suphs0(plat,ng,s_lat%gv,iv)

C --- Loop over sites ---
      iv0 = 0
      allocate(cs(ng))
      allocate(sn(ng))
      do  ib = ib1, ib2
        is = s_site(ib)%spec
        tau = s_site(ib)%pos
        lmxl = s_spec(is)%lmxl
        rg = s_spec(is)%rg
        if (lmxl == -1) cycle
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        if (mode == 2) cofh = 0
        nlm = (lmxl+1)**2
        call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)

        allocate(qkl(nlm*(kmax+1)))
        call pshpr(1)
        call rhogkl(ib,ib,nsp,modgkl,0d0,s_site,s_spec,s_rhat,kmax,qkl)
        call poppr
        call rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,mod(mode/10,10)>=2,qkl,
     .    nlm,ng,g2,yl,cs,sn,cg)
        if (mod(mode,100000) >= 200) call rhgcm3(-z,vol,ng,cs,sn,cg)
        deallocate(qkl)
        iv0 = iv0+nlm
      enddo

      deallocate(yl,g2,iv,cs,sn)
      call tcx('rhgcmp')
      end

      subroutine rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,lcor,qkl,nlm,ng,g2,yl,
     .  cs,sn,cg)
C- Convert G_kL expansion of function centered at a site to PW's
      implicit none
C ... Passed parameters
      integer ng,nlm,kmax
      logical lcor
      double precision ceh,cofh,rfoc,rg,vol,qkl(0:kmax,nlm)
      double precision g2(ng),yl(ng,nlm),cs(ng),sn(ng)
      double complex cg(ng)
C ... Local parameters
      integer i,ilm,l,ll,lmxl,m,k
      double precision aa,cfoc,cvol,gam,gamf,pi,y0,fac,sqkl
      double complex phase,cc

      if (nlm == 0) return
      lmxl = ll(nlm)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      gam = 0.25d0*rg*rg
      gamf = 0.25d0*rfoc*rfoc
      cvol = 1d0/vol
      cfoc = -4d0*pi*y0/vol
      do  i = 1, ng
        phase = dcmplx(cs(i),sn(i))
        aa = dexp(-gam*g2(i))*cvol
        cc = aa*phase*(0d0,1d0)
        ilm = 0
        do  l = 0, lmxl
          cc = cc*(0d0,-1d0)
          do m = -l,l
            ilm = ilm+1
            fac = 1
            sqkl = 0
            do  k = 0, kmax
              sqkl = sqkl + qkl(k,ilm)*fac
              fac = -g2(i)*fac
            enddo
            cg(i) = cg(i) + sqkl*cc*yl(i,ilm)
          enddo
        enddo

        if (lcor) then
          aa = cfoc*dexp(gamf*(ceh-g2(i)))/(ceh-g2(i))
          cg(i) = cg(i) + cofh*aa*phase
        endif

      enddo

      end
      subroutine rhgcm3(z,vol,ng,cs,sn,cg)
C- PW expansion of Z * delta(r)
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :size of delta-function
Ci   vol   :cell volume
Ci   ng    :number of G-vectors
Ci   cs    :cos(-p*G)
Ci   sn    :cos(-p*G)
Co Outputs
Co   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   26 Oct 01
C ----------------------------------------------------------------------

      implicit none
C ... Passed parameters
      integer ng
      double precision z,vol,cs(ng),sn(ng)
      double complex cg(ng)
C ... Local parameters
      integer i
      double complex phase

      do  i = 1, ng
        phase = dcmplx(cs(i),sn(i))
        cg(i) = cg(i) + z*phase/vol
      enddo

      end
