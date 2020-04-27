      subroutine asavqm(mode,s_ctrl,s_pot,s_lat,s_spec,nlmf,vh,qmp,vval)
C- Electrostatic Multipole moments from wave function products
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nclasp nbas nl ips ipc rmax nspin
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ics ipc ips rmax dclabl
Cio    Passed to:  asaqmp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qc qpp pmpol
Cio    Passed to:  asaqmp
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat vol awald nkd nkq nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg jcg indxcg pos qlv dlv symgr ag
Cio    Passed to:  asaqmp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxl lmxa lmxf z
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  1s digit:
Co Outputs
Co   Electrostatic potential
Co   vh
Co   qmp
Co   vval
Cr Remarks
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nlmf
      double precision vval(nlmf,*),qmp(nlmf,*),vh(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      type(str_strxsu) s_strx
      integer, allocatable :: lmxa(:),lmxf(:),lmxl(:)
      real(8), allocatable :: wk(:),yl(:),z(:),pot0l(:),s(:)
C ... Local parameters
      double precision alat,plat(3,3),qlat(3,3),awald,vol,tau(3),xx
      integer ib,ipr,iprint,is,jb,jc,ll,lmax,lmxb,lmxli,nbas,nclasp,
     .  nkd,nkq,nl,nlm0,nlmb,nlmbx,nlml,nlmp,npow,nrx,ic,j
      parameter (nlm0=121)
      double precision dl(2,nlm0),hl(nlm0),rmax,yy,hh,xi(0:20)
      double precision valj(nlm0),pot0(nlm0),q(3)
      character*4 cib

C ... Setup
      nclasp = s_ctrl%nclasp
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      vol = s_lat%vol
      awald = s_lat%awald
      nkd = s_lat%nkd
      nkq = s_lat%nkq
      ipr = iprint()
      if (mode /= 0) call rxi('asavqm: not ready for mode',mode)
C     y0 = 1/sqrt(4*pi)

      allocate(lmxa(nclasp),lmxf(nclasp),lmxl(nclasp))
      allocate(z(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxa',1,lmxa,xx)
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxf',1,lmxf,xx)
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxl',1,lmxl,xx)
      call spec2class(s_spec,nclasp,s_ctrl%ics,'z',1,xx,z)

C     call iinit(lmxl,nclasp)

C --- Multipole moments and cofficients to Hankel expansion ---
      allocate(pot0l(nlmf*nbas))
      call pshpr(ipr-10)
      call asaqmp(0,s_ctrl,s_pot,s_lat,lmxa,lmxf,nlmf,qmp)
      call poppr
      call pvvqm1(nbas,s_ctrl%ipc,lmxf,z,s_pot%qc,nlmf,qmp,pot0l)

C --- For each site, calculate electrostatic at MT boundary ---
C ... Expand about ib from jb
      do  ib = 1, nbas
      write(cib,'(i4)') ib
      is = s_ctrl%ips(ib)
      ic = s_ctrl%ipc(ib)
      lmxli = s_spec(is)%lmxl
      rmax = s_ctrl%rmax(ic)
      lmxli = lmxl(ic)
      nlml = (lmxli+1)**2

C ... Setup for reduced strux
      lmax = 2*ll(nlmf)
      nrx = max(nkd,nkq)
      allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
C ... Setup for strux
      call strxsu(nlml,-99,2*nl-2,1,0,nbas,s_ctrl%ips,s_lat%cg,
     .  s_lat%jcg,s_lat%indxcg,nlmbx,nlmp,npow,s_strx)
      if (nlmp > nlm0) call rxi('asavqm: increase nlm0 to',nlmp)
      allocate(s(nlml*nlmbx))

C ... Loop over all connecting vectors to ib
      do  jb = 1, nbas
C ... lmxb is for now the l-cutoff for the head
C     jc = s_ctrl%ipc(ib)
      jc = s_ctrl%ipc(jb)
      lmxb = lmxf(jc)
      nlmb = (lmxb+1)**2
C ... Connecting vector is pos(jb)-pos(ib) for expansion at ib
      call dpscop(s_lat%pos,tau,3,3*jb-2,1,1d0)
      call dpsadd(tau,s_lat%pos,3,1,3*ib-2,-1d0)
      call shorbz(tau,tau,plat,qlat)
C ... Real reduced strux for this connecting vector
      call dpzero(q,3)
      call hsmqe0(lmxb+lmxli,0d0,110,q,tau,nrx,nlm0,wk,yl,
     .  awald,alat,s_lat%qlv,nkq,s_lat%dlv,nkd,vol,dl)
      call dcopy(nlmp,dl,2,hl,1)
C ... Structure constant matrix
      call hstrux(0d0,nlml,nlmb,nlmp,npow,1,1,s_strx%ikl,
     .  s_strx%jkl(lmxb)%p,s_strx%ip,s_strx%cf,hl,s)
C ... ropbes avoids OKA conventions for Bessel functions
      call ropbes(rmax,0d0,lmxli,yy,hh,xi,1,1)
      call dpscop(pot0l,pot0,nlmb,1+nlmf*(jb-1),1,1d0)
      call dpzero(valj,nlml)
      call sumtl1(nlml,nlmb,s,1,1,pot0,xi,hh,rmax,valj)
      call dpsadd(vval,valj,nlml,1+nlmf*(ib-1),1,1d0)
      enddo
      deallocate(wk,yl,s)
      deallocate(s_strx%cf,s_strx%ip,s_strx%ikl)
      do  j = 0, ll(nlmbx)
        if (associated(s_strx%jkl(j)%p)) then
          deallocate(s_strx%jkl(j)%p)
        endif
      enddo
      enddo
C ... Make vh and Printout
      call pvvqm2(1,nbas,s_ctrl%ipc,lmxl,z,s_ctrl%rmax,
     .  s_pot%qc,nlmf,qmp,pot0l,vval,vh)
      deallocate(lmxa,lmxf,lmxl,z,pot0l)

      end

      subroutine sumtl1(nlml,nlmb,s,nr,nr1,poti,xi,h,rofi,vl)
C- Add bessel tails for potential
      implicit none
      integer nlmb,nlml,nr,nr1
      double precision vl(nr1,nlml),rofi(*),xi(nr,0:*),s(nlml,nlmb),
     .  poti(nlmb),h(*)
C Local
      integer i,ilmb,ilml,l,ll,lmxl,m
      double precision sum

      lmxl = ll(nlml)
      do  i = 1, nr1
        h(i) = 1d0
      enddo
      ilml = 0
      do  l = 0, lmxl
        do  m = -l, l
          ilml = ilml+1
          sum = 0d0
          do  ilmb = 1, nlmb
            sum = sum + s(ilml,ilmb)*poti(ilmb)
          enddo
          do  i = 1, nr1
            vl(i,ilml) = vl(i,ilml) + sum*h(i)*xi(i,l)
          enddo
        enddo
        do  i = 1, nr1
          h(i) = h(i)*rofi(i)
        enddo
      enddo
      end

      subroutine pvvqm1(nbas,ipc,lmxf,z,qc,nlmf,qmp,pot0)
C- Make coffs pot0 from qmom
      implicit none
      integer nbas,nlmf,ipc(*),lmxf(*)
      double precision qmp(nlmf,nbas),pot0(nlmf,nbas),z(*),qc(*)
      integer ib,ilm,k,l,ll,nlm,ic
C     integer ipr,iclbsj,iprint
      double precision df,pi,y0,xx
      character*4 cib

C     ipr = iprint()
      pi = 4*datan(1d0)
      y0 = 1/sqrt(4*pi)

C     if (ipr >= 30) print 649
      do  ib = 1, nbas
        ic = ipc(ib)
        write(cib,'(i4)') ib

        nlm = (lmxf(ic)+1)**2
        do  ilm = 1, nlm
          l = ll(ilm)
          df = 1d0
          do  k = 0, l
            df = df*(2*k+1)
          enddo
          xx = qmp(ilm,ib)
          if (ilm == 1) xx = xx + (qc(ic)-z(ic))*y0
          pot0(ilm,ib) = 2d0*xx*4d0*pi/df
C          if (ipr >= 30 .and. dabs(qmp(ilm,ib)) > 1d-5) then
C            if (iclbsj(ic,ipc,nbas,1) == ib .or. ipr >= 50) then
C            if (ilm == 1) then
C              print 650, cib,ilm,qmp(ilm,ib),z(ic)-qc(ic),pot0(ilm,ib)
C            else
C              print 651, cib,ilm,qmp(ilm,ib),pot0(ilm,ib)
C            endif
C            endif
C          endif
C         cib = ' '
        enddo
      enddo

C  649 format(/'  ib','  ilm',7x,'Qmp',6x,'Z-Qc',7x,'Pot0')
C  650 format(a4,i4,2x,f12.6,f8.3,f12.6)
C  651 format(a4,i4,2x,f12.6,8x,f12.6)

      end

      subroutine pvvqm2(ib1,ib2,ipc,lmxl,z,rmax,qc,nlmf,qmp,pot0,vval,
     .  vh)
C- Printout
      implicit none
      integer ib1,ib2,nlmf,ipc(*),lmxl(*)
      double precision qmp(nlmf,ib2),pot0(nlmf,ib2),vval(nlmf,ib2),
     .  z(*),rmax(*),qc(*),vh(*)
      logical lpr
      integer ib,ilm,ipr,iprint,nlm,ic,iclbsj
      double precision pi,y0,xx,qt
      character*8 cib

      ipr = iprint()
      pi = 4*datan(1d0)
      y0 = 1/sqrt(4*pi)

      if (ipr >= 30) print 1
    1 format('  ib  ic  ilm',7x,'Qmp',6x,'Z-Qc',7x,'Pot0',8x,'Vval',8x,
     .        'Vmad',6x,'Vh(rmt)')
      do  ib = ib1, ib2
        ic = ipc(ib)
        write(cib,'(2i4)') ib,ic

        nlm = (lmxl(ic)+1)**2
        do  ilm = 1, nlm
          xx = qmp(ilm,ib)
          if (ilm == 1) xx = xx + (qc(ic)-z(ic))*y0
          lpr = iclbsj(ic,ipc,ib2,1) == ib .and. ipr >= 30 .or.
     .          ipr >= 50
          if (ilm == 1) then
            qt = qmp(ilm,ib)/y0 + (qc(ic)-z(ic))
            vh(ib) = vval(1,ib)*y0 + 2*qt/rmax(ic)
            if (lpr) print 2,cib,1,qmp(1,ib),z(ic)-qc(ic),pot0(1,ib),
     .        vval(1,ib),vval(1,ib)*y0,vh(ib)
 2          format(a8,i4,2x,f12.6,f8.3,4F12.6)
          else
            if (lpr .and. dabs(qmp(ilm,ib)) > 1d-6) print 3,cib,ilm,
     .          qmp(ilm,ib),pot0(ilm,ib),vval(ilm,ib)
 3          format(a8,i4,2x,f12.6,8x,2F12.6)
          endif
          cib = ' '
        enddo
      enddo

      end
