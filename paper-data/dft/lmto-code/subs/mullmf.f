      subroutine mullmf(nbas,s_site,s_spec,iprmb,z,n,nspc,
     .  iq,isp,mode,nsites,lsites,nchan,lchan,lmdim,nddos,doswt)
C- Make Mulliken decomposition of the norm at this k-point
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxb
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nbas  :size of basis
Ci   z,n   :eigenvectors and dimension n
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   iq    :qp index; used for printout only
Ci   mode  :0 all sites atom-resolved
Ci         :1 all sites l-resolved
Ci         :2 all sites lm-resolved
Ci         :3 site list atom-resolved
Ci         :4 site list l-resolved
Ci         :5 site list lm-resolved
Ci   mode,nsites,lsites,nchan,lchan (see sumlst.f)
Ci   nddos :dimensions doswt
Ci   lmdim :leading dimension of lchan
Co Outputs:
Co   doswt :DOS weights for writing to moms file (lmdos.f)
Cr Remarks
Cr   mullmf sums contributions from orbitals, each into into a
Cr   particular channel identified by lchan.  Each orbital has a site
Cr   index ib and an lm index, which is the information to locates the
Cr   channel index, held by lchan.  lchan contains information for a
Cr   particular (lm,ib), or (l,ib), or (1,ib), depending on mode; thus
Cr   each orbital can be associated with a particular channel, through
Cr   lchan (multiple orbitals may contribute to a single channel).
Cr
Cr   Orthogonal basis : D_in = (z_in)+ z_in where
Cr     i = orbital index and n = band index
Cr   Nonorthogonal basis : D_in = (z~_in)+ z_in
Cr     Here z~ is contravariant form of z.
Cr     Overlap matrix is S = (z z+)^-1
Cr     z~+ = z+ S = z+ (z+)^-1 z^-1 = z^-1
Cu Updates
Cu   26 May 17 doswt no longer initialized
Cu   10 Nov 11 Begin migration to f90 structures
Cu   08 Mar 10 Some bug fixes when number of channels are site dependent
Cu             Altered argument list
Cu   08 Jul 08 Dimension dowst separately from z
Cu   07 Jun 06 Extended to noncollinear case.  Altered argument list
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nspc,iq,isp,n,iprmb(*),mode,nsites,lmdim,
     .  lsites(nsites),nchan,lchan(lmdim,*),nddos
      double precision doswt(nchan,nddos*nspc,nspc)
      double complex z(n,nspc,n*nspc)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical lmp,lp,atp
      integer isite,ib,is,ichan,lmxb,iband,iprint,ipr,iprmin,lgunit,
     .  ispc,ksp,stdo,nev,nx
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      integer io,ikap,l,nlm1,nlm2,ilm,jlm,i,j,norb
      double precision xx
      integer ipiv(n*nspc)
      complex(8),allocatable:: zt(:,:,:),work(:,:,:)
      allocate(zt(n*nspc,n,nspc),work(n*nspc,n,nspc))

C ... Setup
      stdo = lgunit(1)
      call tcn ('mullmf')
C     nx = total hamiltonian dimension
      nx = n*nspc
C     Number of eigenvectors must equal hamiltonian dimemnsion
      nev = nx

C ... Form contravariant (z~)+ = z^-1
      call zcopy(nx**2,z,1,zt,1)
      call zgetrf(nx,nx,zt,nx,ipiv,j)
      if (j /= 0) call rx('mullmf: failed to generate overlap')
      call zgetri(nx,zt,nx,ipiv,work,nx**2,j)

C      call zprm('z',2,z,nx,nx,nx)
C      call zprm('zt',2,zt,nx,nx,nx)

      call sanrg(.true.,mode,0,5,' mullmf:','mode')
      iprmin = 80*iq*isp
      ipr = iprint()
      if (ipr >= iprmin) write(stdo,1)
    1 format (' mullmf:  site spec  lmax  norb  l  nlm1   nlm2  ikappa',
     .        '  offset ilm ichan')
!     call dpzero(doswt,nchan*nddos*nspc)
      atp = .false.
      lp  = .false.
      lmp = .false.
C ... Decompose by atom, l, or lm
      if (mode == 0 .or. mode == 3) atp = .true.
      if (mode == 1 .or. mode == 4) lp  = .true.
      if (mode == 2 .or. mode == 5) lmp = .true.
      if (mode < 3) nsites = nbas

C --- Loop over channels ---
      do  isite = 1, nsites
        ib = lsites(isite)
        is = s_site(ib)%spec
        lmxb = s_spec(is)%lmxb
C   ... Loop over all orbitals centered at this site
        call orbl(ib,0,n,iprmb,norb,ltab,ktab,xx,offl,xx)
C   ... Block into groups of consecutive l
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
        if (ipr >= iprmin) write(stdo,2) ib,is,lmxb,norb
    2   format (9x,i4,1x,i3,6x,i1,3x,i2)

C       In the noncollinear case, isp=1 always => need internal ispc=1..2
C       ksp is the current spin index in both cases:
C       ksp = isp  in the collinear case
C           = ispc in the noncollinear case
C       whereas ispc=1 for independent spins, and spin index when nspc=2
        do  ispc = 1, nspc
        ksp = max(ispc,isp)

        do  io = 1, norb
          l  = ltab(io)
          ikap  = ktab(io)
          nlm1 = l**2+1
          nlm2 = (l+1)**2
C         i = orbital index in iprmb order
          i = offl(io)
          if (i > n) call rx (' bug in mullmf: i>n')
          if (ipr >= iprmin) write(stdo,3) l,nlm1,nlm2,ikap,i+1
    3     format (33x,i1,3x,i2,5x,i2,6x,i1,3x,i5)
          do  ilm = nlm1, nlm2
            i = i + 1
            if (atp) jlm = 1
            if (lp)  jlm = l+1
            if (lmp) jlm = ilm
            if (jlm > lmdim) cycle
            ichan = lchan(jlm,isite)
C            call mchan(lmdim,0d0,0d0,0d0,0d0,0,nsites,0,isite,jlm,1,
C     .        ichan,lchan)
            if (ichan > nchan) call rx(' bug in mullmf: ichan>nchan')
            if (ichan == 0) cycle
            if (ipr >= iprmin) write(stdo,4) ilm,ichan
    4       format (65x,i2,1x,i4)
            do  iband = 1, nev
              doswt(ichan,iband,ispc) = doswt(ichan,iband,ispc) +
     .          dble( zt(iband,i,ispc)*z(i,ispc,iband) )
            enddo
          enddo
        enddo
        enddo

      enddo

C ... debugging ... xx should be 1
C      do  iband = 1, nev
C        xx = 0
C        do  ispc = 1, nspc
C        do  ichan = 1, nchan
C          xx = xx + doswt(ichan,iband,ispc)
C        enddo
C        enddo
C        print *, iband, sngl(xx)
C      enddo

      deallocate(zt,work)
      call tcx('mullmf')

      end
