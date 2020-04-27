      subroutine grp2av(s_ctrl,s_site,s_spec,s_rhat)
C- Average l=0 parts of site densities that share common GRP2
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: nclass ics ipc
Co     Stored:    *
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Cio    Passed to: pgrp2a prrhat
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: nr z a lmxl grp2 rmt
Co     Stored:    *
Cio    Passed to: pgrp2a prrhat
Cio  s_rhat
Ci     Elts read: rho1 rho2
Co     Stored:    *
Co     Stored:    rho1 rho2 averaged on output according to common GRP2
Cio    Passed to: pgrp2a prrhat
Ci Inputs
Co Outputs
Cl Local variables
Cl   jgroup:jgroup(i) = abs(igroup(i)), i<nclass
Cl         :jgroup(nclass+i) = class order:
Cl         :classes are permuted by increasing group index.  Thus:
Cl         :1 classes associated with the same group index are contiguous
Cl         :2 Actual class for ith element is jgroup(nclass+i1)
Cl         :  using convention ic = 0:nclass-1
Cr Remarks
Cr   Site densities that share a common common nonzero value of
Cr   s_spec%grp2 are averaged.  More precisely the l=0 part of rho1 and
Cr   rho2 are averaged (with possible spin flips on some sites, as
Cr   described below).  Only the l=0 part is taken because it is
Cr   independent of rotations.
Cr
Cr   The operation is intended to enforce symmetry between sites that
Cr   might be equivalent, or nearly so but are not equivalent according
Cr   to the given symmetry operations.  For example a 2-atom supercell
Cr   cell of the bcc structure, are equivalent but not required to be so
Cr   according to the 48 point group operations.  Or, magnetic atoms
Cr   of antiferromagnet could be equivalent apart from a spin flip.
Cr
Cr   For the operation to be meaningful, the atomic number Z and radial
Cr   mesh should be common to all sites averaged.  This routine will
Cr   stop if sites with different Z or number of radial mesh points are
Cr   different.
Cr
Cr   All sites with a common GRP2 (specified by s_spec%grp2
Cr   corresponding to the site) are combined.  In addition sites with a
Cr   common -GRP2 are included in the average with spin flips.  The l=0
Cr   part of rho1 and rho2 are averaged, and then redistributed over
Cr   sites.
Cu Updates
Cu   31 Jul 13 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: igroup(:),jgroup(:)
      integer, pointer :: ips(:)
      real(8), allocatable:: rho1(:,:),rho2(:,:),qs(:,:),rofi(:),rwgt(:)
C ... Local parameters
      integer cgrp,grps,i,i0,ib,intopt,is,isp,ispc,k,
     .  nbas,ngr2,nlml,nr,nsp,stdo
      integer nglob,iprint,iasum
      double precision a,ddot,rmt,srfpi,xx,z

C --- Make array of grp2, one element for each class ---
      nbas = nglob('nbas')
      allocate(igroup(0:nbas-1))
      ips => s_ctrl%ips
      call spec2class(s_spec,nbas,ips,'grp2',1,igroup,xx)
      if (iasum(nbas,igroup,1) == 0) then
        deallocate(igroup)
        return
      endif

      call tcn('grp2av')
      call info0(30,1,0,
     .  ' Average l=0 site densities with common GRP2 ...')
      stdo = nglob('stdo')
      intopt = 10*nglob('lrquad')
      nsp = nglob('nsp')
      srfpi = dsqrt(16d0*datan(1d0))

C --- Copy iabs(igroup) to jgroup
      allocate(jgroup(0:2*nbas-1),qs(nbas,2))
      call dpzero(qs,nbas*2)
      do  ib = 0, nbas-1
        jgroup(ib) = iabs(igroup(ib))
      enddo

C --- Sort groups and find first site with grp2>0 ---
      call ivshel(1,nbas,jgroup,jgroup(nbas),.true.)
      do  ib = 0, nbas-1
        i0 = ib
        if (jgroup(jgroup(nbas+ib)) /= 0) exit
      enddo
      if (i0 == nbas-1) call rx('bug in grp2av')

C --- Average rho by grp2 ---
      ngr2 = 0  ! Number of groups
      grps = 0  ! Number of sites in current group
      cgrp = 0  ! Group number
C ... Loop through classes in order of permutation table
C     Loop order ensures sites with common group are contiguous in table.
      do  i = i0, nbas-1
C   ... New group if group ID changes
        if (cgrp /= jgroup(jgroup(nbas+i))) then
          if (grps /= 0)      ! Distribute prior grp2 before starting new one
     .      call pgrp2a(s_site,s_spec,i,ngr2,grps,nbas,nr,nsp,igroup,
     .      jgroup,rho1,rho2,qs,s_rhat)

          cgrp = jgroup(jgroup(nbas+i))
          ib = jgroup(nbas+i)+1 ! Site for first member of group
          is = ips(ib)
          grps = 0
          ngr2 = ngr2+1
          nr = s_spec(is)%nr
          z = s_spec(is)%z
          if (allocated(rho1)) deallocate(rho1,rho2,rofi,rwgt)
          allocate(rho1(nr,nsp),rho2(nr,nsp),rofi(nr),rwgt(nr))
          call dpzero(rho1,nr*nsp)
          call dpzero(rho2,nr*nsp)
          rmt = s_spec(is)%rmt; a = s_spec(is)%a
          call radmsh(rmt,a,nr,rofi)
          call radwgt(intopt,rmt,a,nr,rwgt)
        endif
        ib = jgroup(nbas+i)+1 ! True site index for this site
        is = s_site(ib)%spec
C   ... Sum rho1,rho2 of this site into group average
        if (nr /= s_spec(is)%nr .or. z /= s_spec(is)%z) call
     .    rx('grp2av: illegal group association: mismatch in nr or z')
        do  isp = 1, nsp
          ispc = isp
          if (igroup(ib-1) < 0) ispc = nsp+1-isp
          nlml = (s_spec(is)%lmxl+1)**2; k = nlml*(isp-1)
          call daxpy(nr,1d0,s_rhat(ib)%rho1(1,1+k),1,rho1(1,ispc),1)
          call daxpy(nr,1d0,s_rhat(ib)%rho2(1,1+k),1,rho2(1,ispc),1)
          qs(ib,1) = qs(ib,1) +
     .      srfpi*ddot(nr,s_rhat(ib)%rho1(1,1+k),1,rwgt,1)
          qs(ib,2) = qs(ib,2) +
     .      srfpi*ddot(nr,s_rhat(ib)%rho2(1,1+k),1,rwgt,1)
        enddo
        grps = grps+1
      enddo

C ... Distribute last grp2
      call pgrp2a(s_site,s_spec,nbas,ngr2,grps,nbas,nr,nsp,igroup,
     .  jgroup,rho1,rho2,qs,s_rhat)
      if (allocated(rho1)) deallocate(rho1,rho2,rofi,rwgt)
      deallocate(igroup,jgroup,qs)
      if (iprint() > 50) call prrhat(nbas,s_site,s_spec,s_rhat)
      call tcx('grp2av')

      end

      subroutine pgrp2a(s_site,s_spec,i,ngr2,grps,nbas,nr,nsp,
     .  igroup,jgroup,rho1,rho2,qs,s_rhat)
C- Redistribute avg rho(l=0) into each site associated with current group
C ----------------------------------------------------------------------
Ci Inputs
Ci   i     :Counter marking site after last site in current group
Ci         :True site ib is given by jgroup(nbas+i)+1; see jgroup below
Ci         :Note i ranges from 0..nbas-1, not 1,nbas
Ci   ngr2  :Number of groups
Ci   grps  :Number of sites in group
Ci   nbas  :size of basis
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   igroup:group associated with each class, derived from s_spec%grp2
Ci   jgroup:jgroup(i) = abs(igroup(i)), i<nbas
Ci         :jgroup(nbas+i) = site order:
Ci         :sitees are permuted by increasing group index.  Thus:
Ci         :1 sitees associated with the same group index are contiguous
Ci         :2 Actual site for ith element is ib = jgroup(nsite+i) + 1
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Cio    Passed to: *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxl
Co     Stored:    *
Cio    Passed to: *
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Ci   qs    :sphere charges by site before averaging (for printout only)
Cio  s_rhat
Ci     Elts read: rho1 rho2
Co     Stored:    *
Cio    Passed to: *
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   31 Jul 13
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer i,ngr2,grps,nbas,nr,nsp
      integer igroup(0:nbas-1),jgroup(0:2*nbas-1)
      double precision rho1(nr,nsp),rho2(nr,nsp),qs(nbas,2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(nbas)
      type(str_spec)::  s_spec(*)
      type(str_rhat)::  s_rhat(nbas)
C ... Local parameters
      integer i1,iprint,isp,ispc,stdo,nglob,ib,is,k,nlml
      double precision qbar,q2,rms

      stdo = nglob('stdo')
      call dscal(nr*nsp,1/dble(grps),rho1,1)
      call dscal(nr*nsp,1/dble(grps),rho2,1)

C --- Redistribute rho1,rho2 of this group into each class ---
      qbar = 0; q2 = 0
      do  i1 = i-grps, i-1
        ib = jgroup(nbas+i1)+1 ! True site index
        is = s_site(ib)%spec
        nlml = (s_spec(is)%lmxl+1)**2
        do  isp = 1, nsp
          ispc = isp
          if (igroup(ib-1) < 0) ispc = nsp+1-isp
          k = nlml*(isp-1)
          call dcopy(nr,rho1(1,ispc),1,s_rhat(ib)%rho1(1,1+k),1)
          call dcopy(nr,rho2(1,ispc),1,s_rhat(ib)%rho2(1,1+k),1)
        enddo
        qbar = qbar + qs(ib,1)
        q2 = q2 + qs(ib,1)**2
      enddo
      qbar = qbar/grps; q2 = q2/grps; rms = dsqrt(dmax1(q2-qbar**2,0d0))
      if (iprint() >= 30) write(stdo,2) grps,ngr2,qbar,rms
    2 format(i5,' elts in group',i3,'  <Q>=',f9.6,'  RMS GRP2 DQ=',f9.6)
      end
