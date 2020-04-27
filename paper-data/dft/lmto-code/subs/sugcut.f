      subroutine sugcut(mode,nspec,s_spec,alat,ng,gv,tol)
C- Find max recip for each spec and orbital block, store in struct.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name ngcut lmxa lmxb pz orbp
Co     Stored:     ngcut orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Ci Inputs
Ci   mode  :1 make cutoffs for standard envelope functions
Ci         :2 make cutoffs for extended local orbitals
Ci         :3 combination 1+2
Ci         :10s digit
Ci         :1 return gtop = global gmax that meets tol criterion;
Ci            returns gtop in gv(1,1)
Ci   nspec :number of species
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors (not used for mode>=10)
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   tol   :tolerance in wave function precision:
Ci         :Here wave function cutoff for g-vectors whose wave
Ci          functions are less than tol.
Co Outputs
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 May 09 New 10s digit mode
Cu   16 Aug 04 New mode for getting cutoffs, local orbs.
Cu             Changed argument list
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu    9 May 00 Adapted from nfp su_gvcut.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nspec,ng
      double precision alat,gv(ng,3),tol
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer nkap0,n0
      parameter (nkap0=4,n0=10)
      integer ncut(n0,nkap0),lh(nkap0)
      integer stdo,lgunit,ipr,iprint,is,irep,icut,i,ik,l,lcut,nkapi,
     .  nkap1,nkap2,mode0,nglob
      double precision rsmh(n0,nkap0),eh(n0,nkap0),tpi,tpiba2,gg0,gg,
     .  e,rsm,gam,gmax,top,gtop
      character*8 spid
      character*1 ccc,ccl

      if (ng == 0 .and. mode < 10) return

C     If no extended local orbitals, exit
      if (mode == 2) then
        i = nglob('lpz')
        if (mod(i/10,10) == 0) return
      endif

      stdo = lgunit(1)
      ipr = iprint()
      tpi = 8d0*datan(1d0)
      tpiba2 = (tpi/alat)**2
      mode0 = mod(mode,10)
      if (ipr >= 20) then
C       if (mode == 1) write(stdo,887) tol
C       if (mode == 2) write(stdo,888) tol
        call info2(20,1,0,' sugcut:  make '//
     .    '%?#n>9#global maximum #'//
     .    'orbital-dependent #%-1j'//
     .    '%?#n==1#reciprocal vector cutoffs for##%-1j'//
     .    '%?#n==2#cutoffs for local orbitals,##%-1j'//
     .    '%?#n==3#cutoffs for valence+local orbitals,##%-1j'//
     .    '%?#n==11#reciprocal vector cutoff for##%-1j'//
     .    '%?#n==12#cutoff for local orbitals,##%-1j'//
     .    '%?#n==13#cutoff for valence+local orbitals,##%-1j'//
     .    '%j tol=%,2;3g',mode,tol)
        write(stdo,774)
C  887   format(/' sugcut:  make orbital-dependent reciprocal vector',
C     .    ' cutoffs for tol=',1p,e9.2)
C  888   format(/' sugcut:  orbital-dependent cutoffs for local',
C     .    ' orbitals, tol=',1p,e9.2)
      endif
      call sanrg(.true.,mode0,1,2,'sugcut','1s digit mode')

      gtop = -1
      gg = -1
      do  is = 1, nspec
C       call spacks(0,'spec name',sspec,spid,is,is)
        spid = s_spec(is)%name
        nkap1 = 1
        call uspecb(0,1,s_spec,is,is,lh,rsmh,eh,nkapi)
        if (mode0 > 1) then
          call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkap2)
          if (mode0 == 2) nkap1 = nkapi+1
        else
          nkap2 = nkapi
        endif
        call iinit(ncut,n0*nkap0)
        if (mode0 == 2) ncut = s_spec(is)%ngcut
        gg0 = gg
        do  ik = nkap1, nkap2
        lcut = -1
        do  l  = 0, lh(ik)
          e = eh(l+1,ik)
          rsm = rsmh(l+1,ik)
          if (rsm /= 0) then
          if (l < lh(ik) .and. l > lcut) then
            lcut = l-1
   12       lcut = lcut+1
            if (lcut < lh(ik)) then
            if (rsmh(lcut+2,ik) == rsm .and. eh(lcut+2,ik) == e) goto 12
            endif
          endif
C     ... Get cutoff radius where exp(-gam*gmax)*gmax**l equals tol
          gam = rsm*rsm/4d0
          gmax = 1d0
          do  irep = 1, 10
            gmax = dsqrt(-dlog(tol/gmax**l)/gam)
C           write(stdo,895) irep,gmax
C 895       format('irep,gmax=',i5,f12.6)
          enddo
          if (mode >= 10) then
            gtop = max(gtop,gmax)
          else
C     ... Find first longer vector, icut is one less
          icut = 1
          do  i = 1, ng
            gg = tpiba2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
            if (gg > gmax*gmax) goto 90
            icut = i
            gg0 = gg
          enddo
   90     continue
          top = dexp(-gam*gg0)*dsqrt(gg0)**l
          ccc = ' '
          if (icut == ng) ccc = '*'
          ccl = ' '
          if (l < lcut) ccl = '*'
          if (ipr >= 20) write(stdo,773) spid,l,ccl,rsm,e,gmax,
     .      top,icut,ccc
  773     format(2x,a,i2,a1,f7.2,f7.2,f8.3,1p,e12.2,0p,i8,a)
  774     format(' spec      l    rsm',4x,'eh',5x,'gmax',
     .       4x,'last term',3x,'cutoff')
          ncut(l+1,ik) = icut
          endif
          endif
        enddo
        enddo
        s_spec(is)%ngcut = ncut
      enddo

      if (mode >= 10) then
        gv(1,1) = gtop
      endif

      call info0(30,0,0,' ')

      end
