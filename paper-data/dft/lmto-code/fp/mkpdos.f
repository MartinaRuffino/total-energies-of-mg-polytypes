      subroutine mkpdos(mode,s_site,s_spec,isp,nsp,nspc,nlmax,
     .  ndham,nphimx,nev,nchan,lchan,lmdim,npd,lsites,nsites,ppnl,aus,doswt)
C- Decomposition of norm from projection of w.f into MT sphere
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa rmt
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 or 3 DOS resolved by site R
Ci         :1 or 4 DOS resolved by R and l
Ci         :2 or 5 DOS resolved by R,l,m
Ci         :10s digit (not implemented)
Ci         :1 doswt = separate phi,phidot contributions to norm, i.e.
Ci         :          <phi | evec> ;  <phidot | evec>
Ci         :          In this case, ndp should be at least 2.
Ci         :100s digit
Ci         :0  aus = (phi,phidot) of w.f. at MT boundary
Ci         :1  aus = (val,slo) of w.f. at MT boundary (not implemented)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus : must be at least hamiltonian rank
Ci   nphimx:dmensions aus: global max # of partial waves of a given l
Ci   nev   :number of eigenvectors to accumulate DOS
Ci   nchan :total number of channels
Ci   npd   :second dimension of doswt.
Ci   mode,nsites,lsites,nchan,lchan (see sumlst.f)
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci   aus   :values of (phi,phidot) MT sphere boundary; see makusq
Co Outputs
Co   doswt :dos weights for this qp are added to doswt
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   07 Mar 10 Redesigned to restrict pdos to channels specified by lchan
Cu   13 Feb 02 Extended to local orbitals.
Cu   28 Mar 01 written by MvS.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isp,nsp,nspc,nlmax,ndham,nphimx,nchan,nsites,lsites(nsites),nev,
     .  npd,lmdim,lchan(lmdim,*)
      integer n0,nppn
      parameter (n0=10,nppn=12)
      double precision ppnl(nppn,n0,nsp,*),doswt(nchan,npd,nev,nspc)
      double complex aus(nlmax,ndham*nspc,nphimx,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer lmxa,lmxax,ichan,isite,ib,is,iv,ilm,l,m,ll,jlm,ispc,ksp
      double precision s11,s22,s33,s12,s13,s23,sum,rmt
      double complex au,as,az

      call sanrg(.true.,mode,0,5,' mkpdos','mode')
      call sanrg(.true.,nev,0,ndham*nspc,' mkpdos','nev')
      lmxax = ll(nlmax)

C     ichan = 0
      do  isite = 1, nsites
        if (lsites(1) <= 0) then
          ib = isite
        else
          ib = lsites(isite)
        endif
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        rmt = s_spec(is)%rmt
        lmxa = min(lmxa,lmxax)
        if (lmxa == -1) cycle

C       call pp2hvs(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)
C       call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)

        ilm = 0
        do  l = 0, lmxa
        do  m = -l, l
          ilm = ilm+1
          if (mode == 0 .or. mode == 3) jlm = 1
          if (mode == 1 .or. mode == 4) jlm = l+1
          if (mode == 2 .or. mode == 5) jlm = ilm
          if (jlm > lmdim) cycle
          ichan = lchan(jlm,isite)
          if (ichan == 0) cycle
          if (ichan > nchan) call rx(' bug in mkpdos: ichan>nchan')

          do  ispc = 1, nspc
          ksp = max(ispc,isp)
          az = 0
          do  iv = 1, nev

            au = aus(ilm,iv,1,ksp,ib)
            as = aus(ilm,iv,2,ksp,ib)
            if (nphimx > 2) az = aus(ilm,iv,3,ksp,ib)

C           If (au,as) are coefficients to (u,s), use this
C           s11 = dconjg(au)*au*sab(1,l+1,ksp)
C           s12 = 2*dconjg(au)*as*sab(2,l+1,ksp)
C           s22 = dconjg(as)*as*sab(4,l+1,ksp)

C           If (au,as) are coefficients to (phi,phidot), use this
            s11 = dconjg(au)*au*ppnl(2,l+1,ksp,ib)
            s22 = dconjg(as)*as*ppnl(7,l+1,ksp,ib)
            s33 = dconjg(az)*az*ppnl(8,l+1,ksp,ib)
            s12 = 0
            s13 = 2*dconjg(au)*az*ppnl(9,l+1,ksp,ib)
            s23 = 2*dconjg(as)*az*ppnl(10,l+1,ksp,ib)

            sum = s11+s22+s33 + s12+s13+s23

            doswt(ichan,1,iv,ispc) = doswt(ichan,1,iv,ispc) + sum

          enddo
          enddo
        enddo
        enddo
      enddo

      end
