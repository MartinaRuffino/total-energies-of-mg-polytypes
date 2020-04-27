      subroutine asaqmp(iopt,s_ctrl,s_pot,s_lat,lmxa,lmxf,nlmf,qmp)
C- Multipole moments from 2nd gen ASA wave function products
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nclasp nbas nl nspin
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ipc dclabl
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qpp pmpol
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:jcg indxcg cg pos symgr ag
Cio    Passed to:  *
Ci Inputs
Ci   iopt  :not used: must be zero
Ci   lmxa  :array of augmentation l-cutoffs
Ci   lmxf  :array of l-cutoffs for qpp
Ci   nlmf  :leading dimension and global L-cutoff for qmp
Co Outputs
Ci   qmp   :multipole moments
Cr Remarks
Cr   For nonspherical multipole moments, let
Cr        I^(pp)_l'l''m = int (phi_l' phi_l'' r^m)
Cr        I^(pd)_l'l''m = int (phi_l' phidot_l'' r^m)
Cr        I^(dd)_l'l''m = int (phidot_l' phidot_l'' r^m)
Cr   Then the multipole moments inside sphere R are
Cr     q_M = sum_L',L'' CG_ML'L''
Cr           [ (1+oh)z+_RL' I^(pp)_l'l''m (1+oh)z_RL''   +
Cr             (1+oh)z+_RL' I^(pd)_l'l''m (hz)_RL'' + h.c. +
Cr                (hz)+_RL' I^(dd)_l'l''m (hz)_RL'' ]
Cr         = sum_L',L'' CG_ML'L''
Cr           [ q_pp(L',L'') I^(pp)_l'l''m +
Cr           2*q_pd(L',L'') I^(pd)_l'l''m
Cr             q_dd(L',L'') I^(dd)_l'l''m ]
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   08 Nov 07 (J. Xu) qpp is complex
Cu   23 Aug 01 adapted from makqmp.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer iopt,nlmf,lmxa(*),lmxf(*)
      double precision qmp(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      integer, allocatable :: ipa(:)
      real(8), allocatable :: posc(:),qwk(:),sym(:)
C ... Local parameters
      double precision alat,plat(3,3),qlat(3,3)
      integer ngrp,nrclas,nclasp,ic,nlml,nl,nn,ipr,iprint,
     .  nlmx,lgunit,i,nbas,nqpp,nsp
      character*8 clabl

C ... Setup
      nclasp = s_ctrl%nclasp
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      ngrp = s_lat%nsgrp
      ipr = iprint()
      if (iopt /= 0) call rx('asaqmp: bad  iopt')

C --- Multipole moments at all sites from the qpp ---
      i = nl**2
      nqpp = (i*(i+1))/2
      call qpp2mp(nqpp,nl,nsp,nbas,nlmf,s_ctrl%ipc,lmxa,
     .  s_lat%jcg,s_lat%indxcg,s_lat%cg,s_pot%qpp,s_pot%pmpol,qmp)
c     call psymqp(nlmf,nlmf,1,nbas,qmp)

      if (ipr >= 20)
     .  call awrit1('%N ASAQMP: Make and symmetrize multipole moments'
     .  //' for %i classes',' ',80,lgunit(1),nclasp)

      allocate(ipa(nbas),posc(3*nbas),qwk(nlmf))

C --- For each class, do ---
      do  ic = 1, nclasp
        nlml = (lmxf(ic)+1)**2
        nlmx = max0(nlml,4)

C   ... Make nrclas,ipa,posc
        call psymr0(-2,ic,nbas,s_ctrl%ipc,s_lat%pos,posc,ipa,nrclas)
        do  i  = 1, nrclas
          ipa(i) = (ipa(i)-1)*nlmf
        enddo

C   ... Symmetrize qmp for members of this class
        allocate(sym(nlmx*nlmx*nrclas))
        call symqmp(nrclas,nlml,nlmx,plat,posc,ngrp,s_lat%symgr,
     .    s_lat%ag,qwk,ipa,sym,qmp,nn)
        deallocate(sym)

        call r8tos8(s_ctrl%dclabl(ic),clabl)
        if (ipr >= 30)
     .    call awrit3('%N Class '//clabl//'%a: %16p'//
     .    'nrc = %i,  %i nonspherical elements of %i ',' ',80,lgunit(1),
     .    nrclas,nn,nlml)

C        if (ipr >= 30 .and. ipr < 50)
C     .    call psymqp(nlmf,nlmf,ipa,ipa,qmp)

      enddo

      deallocate(ipa,posc,qwk)

      if (ipr >= 50) call psymqp(nlmf,nlmf,1,nbas,qmp)

      end

      subroutine psymqp(nlmf,nlml,ib1,ib2,qmp)
C- Printout
      implicit none
      integer nlmf,nlml,ib1,ib2
      double precision qmp(nlmf,ib2),fpi,y0
      integer j1,ib,ilm

      print 1
    1 format ('  ib   ilm      qmom',8x,'Qval')
      fpi = 16d0*datan(1d0)
      y0 = 1d0/dsqrt(fpi)
      j1 = 1
      do  ib = ib1, ib2
        print 2,ib,1,qmp(j1,ib),qmp(j1,ib)/y0
    2   format(i4,i6,f12.6,f12.6,2F9.2)
        do  ilm = 2, nlml
          if (dabs(qmp(ilm,ib)) > 1d-6) print 3,ilm,qmp(ilm,ib)
    3     format(4x,i6,f12.6)
        enddo
      enddo

      end


      subroutine qpp2mp(nqpp,nl,nsp,nbas,nlmf,ipc,lmxa,
     .  jcg,indxcg,cg,qpp,pmpol,qmp)
C- Make multipole moments from the qpp and pmpol
      implicit none
      integer nqpp,nbas,nsp,nl,nlmf,lmxa(*),ipc(nbas),jcg(*),indxcg(*)
      double precision qmp(nlmf,nbas),pmpol(nl,nl,2*nl-1,3,nsp,1),cg(1)
      double complex qpp(nqpp,4,nsp,nbas)
C Local
      integer ib,nlm,iqpp,iqpd,ilm1,l1p1,ll,ilm2,l2p1,ix,icg,mlm,lm,isp,ic
      double precision wij

      call dpzero(qmp,nlmf*nbas)

C --- For each site, make qmp ---
      do  isp = 1, nsp
      do  ib = 1, nbas
        ic = ipc(ib)
        nlm = (lmxa(ic)+1)**2

C   ... debugging
C        iqpp = 0
C        snot = 0
C        do  ilm1 = 1, nlm
C        do  ilm2 = 1, nlm
C          if (ilm2 <= ilm1) then
C            wij = 2
C            if (ilm1 == ilm2) wij = 1
C            iqpp = iqpp+1
C            print *, ilm1,ilm2,iqpp,qpp(iqpp,1,isp,ib)
C            snot(ilm1,ilm2) = qpp(iqpp,1,isp,ib)
C          endif
C        enddo
C        enddo
C        print *, 'ib,isp=',ib,isp
C        call prmx('qpp',snot,16,16,16)

        iqpp = 0
        iqpd = 0
        do  ilm1 = 1, nlm
        l1p1 = ll(ilm1)+1
        do  ilm2 = 1, nlm
        l2p1 = ll(ilm2)+1

          ix = max0(ilm1,ilm2)
          ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)

C     ... phi-phi and dot-dot terms
          if (ilm2 <= ilm1) then
            wij = 2
            if (ilm1 == ilm2) wij = 1
            iqpp = iqpp+1
            do  icg = indxcg(ix), indxcg(ix+1)-1
              mlm = jcg(icg)
              if (mlm <= nlmf) then
                lm = ll(mlm)
                qmp(mlm,ib) = qmp(mlm,ib) + cg(icg)*wij*
     .      (dble(qpp(iqpp,1,isp,ib))*pmpol(l1p1,l2p1,lm+1,1,isp,ic) +
     .       dble(qpp(iqpp,4,isp,ib))*pmpol(l1p1,l2p1,lm+1,3,isp,ic))
              endif
                enddo
          endif

C     ... phi-dot terms
          iqpd = iqpd+1
          if (iqpd > nqpp) cycle  !! May be a bug ... needs revisiting
          do  icg = indxcg(ix), indxcg(ix+1)-1
            mlm = jcg(icg)
            if (mlm <= nlmf) then
              lm = ll(mlm)+1
C              if (mlm == 2 .and. ib == 1) then
C                call awrit8(
C     .            'l1,m1= %i %i l2,m2= %i %i iqpd=%i %32p%g %g %g',
C     .            ' ',80,6,l1p1-1,ilm1,l2p1-1,ilm2,iqpd,
C     .            qpp(iqpd,2,isp,ib),cg(icg),
C     .            pmpol(l1p1,l2p1,lm,2,isp,ic))
C              endif
              qmp(mlm,ib) = qmp(mlm,ib) + cg(icg)*2*
     .    (dble(qpp(iqpd,2,isp,ib))*pmpol(l1p1,l2p1,lm,2,isp,ic))

            endif
          enddo
          enddo
          enddo

        enddo
      enddo

      end
