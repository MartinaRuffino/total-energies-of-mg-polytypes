      subroutine smvxt(s_site,s_spec,s_lat,s_pot,nbas,k1,k2,k3,qmom,
     .  gpot0,hpot0,smrho,smpot,f,rhvxt,rcvxt,sgpe,eps0,ieps0)
C- Interaction of an external potential with the smooth density.
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
Ci   nbas  :size of basis
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   smrho :smooth density on real-space mesh, without compensating gaussians
Co Outputs (see also Remarks)
Co   gpot0 :integrals of (compensating gaussians g_RL) * phiext
Co         :is added to the original gpot0 made in vesgcm
Co         :NB: There is a local analog of gpot0 generated in locpt2.
Co         :    which should cancel this term
Co   hpot0 :integrals of semicore smooth Hankels * phiext
Co   smpot :smooth potential phiext is added to smpot
Co   f     :electrostatic contribution to force.
Co   rhvxt :integral n0~ phiext
Co         :where n0~ is sm. density including compensating gaussians
Co   sgpe  :sum_RL integral qmom_RL g_RL phiext = integral [n0~-n0] phiext
Co   rcvxt :integral (qcorg-z + rhoc) phiext
Cl Local variables
Cl   smq   :integral of smooth density n0
Cl   gpote :integrals of compensating gaussians g_RL * phiext
Cl   u00   :integral n0 (phi[n0]+phiext) = n0 (phi0+phiext)
Cl   u0g   :integral n0 [phi0~-phi0)
Cl   ugg   :integral [n0~-n0] [phi0~-phi0]
Cl         :ugg is not used
Cs Command-line switches
Cs   --wsmpot : Write smooth potential to file
Cr Remarks
Cr  This routine is an analog of smves.
Cr  Here the potential is an external phiext; it is not generated from the density
Cr  This routine:
Cr    * adds phiext to smpot
Cr    * Add integrals (compensating gaussians g_RL) *  phiext to gpot0
Cr    * Add integrals (semicore smooth Hankels) *  phiext to hpot0
Cr    * Counts contribution to the force
Cr    * Computes terms for total energy
Cr   Let n0  = smooth density without the compensating sum_RL Q_RL g_RL
Cr       n0~  = n0 + sum_RL Q_RL g_RL
Cr      g_RL  = gaussian in RL channel
Cr      h_R   = l=0 sm hankel in RL channel, (represents core densities)
Cr    qmom_RL = multipole moment in RL channel of (n_R(r) - n0_R(r))
Cr              so that int n_RL(r)-n0_RL(r) = qmom_RL * g_RL(r)
Cr      gpot0 = vector of integrals g_RL * [phi0~ + phiext]
Cr              The latter term is added to gpot0 [vesgcm]
Cr      hpot0 = integral h_R * (phi0~+phiext) (contributions from core)
Cr              The latter term is added to hpot0 [vesgcm]
Cr       ue0  = integral n0 phiext
Cr       sgpe = integral [n0~-n0] phiext
Cr     rhovxt = ue0 + sgpe  [valence contribution to rho*phiext]
Cr     rhcvxt = integral (core density) * phiext  [corresponding core contribution]
Cu Updates
Cu   16 Jun 16 Adapted from smves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer k1,k2,k3,nbas,ieps0
      double precision rhvxt,sgpe,rcvxt,eps0
      double precision qmom(*),f(3,nbas),gpot0(*),hpot0(nbas)
      double complex smrho(k1,k2,k3,2),smpot(k1,k2,k3,2)
C ... For structures
!       include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
C ... Dynamically allocated local arrays
      real(8), allocatable :: gpote(:),hpote(:)
      complex(8), allocatable :: cv(:,:),cgsum(:)
C     complex(8), allocatable :: cn(:)
C ... Local parameters
      integer ib,ilm,ipr,is,iv0,lfoc,lgunit,lmxl,n1,n2,
     .  n3,ng,ngabc(3),nglob,nlm,nsp,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision ceh,cofg,cofh,pi,qcorg,qcorh,qsc,
     .  rfoc,srfpi,ue0,vol,y0,z,qsmc,s1,s2,zsum
      procedure(integer) :: iprint
C     procedure(complex(8)) :: zdotc

C ... Setup
      call tcn('smvxt')
      ipr   = iprint()
      stdo  = lgunit(1)
      nsp   = nglob('nsp')
      pi    = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      y0    = 1d0/srfpi
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol

      call info0(30,1,-1,' smvxt: add external potential')
      if (nsp == 2) call rx('smvxt not ready for nsp=2')

C ... Contributions from vext to smpot, gpot0, hpot0, total energy
      allocate(cv(ng,nsp)) ! For now, cv can only be one spin.
      call dcopy(2*ng*nsp,s_pot%smvextc,1,cv,1)

C ... Smooth density, potential in reciprocal space; add phiext to smpot
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,-1)
      call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,-1)

C ... add phiext to smpot
C     allocate(cn(ng))
C     call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smrho,cn)
C     ue0 = vol*zdotc(ng,cn,1,cv,1)    ! integral smrho * vext (sans compensating gaussians)
      call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smrho,smpot,s1,s2)  ! smrho * ves
      call gvaddf(ng,s_lat%kv,k1,k2,k3,cv,smpot)
      call mshdot(vol,1,n1,n2,n3,k1,k2,k3,smrho,smpot,ue0,s2) ! smrho * (ves+phiext)
      ue0 = (ue0-s1)*n1*n2*n3          ! integral smrho * vext (smrho sans compensating gaussians)

C ... Analysis of eps
      if (ieps0 > 0) then
        call info0(15,0,0,' ... analysis of screening')
        call smvxte(10+ieps0,s_lat,s_pot,n1,n2,n3,eps0,smpot,smrho)
      else
        call info0(30,0,0,' ... ')
      endif

C#ifdefC DBGVXT
C      call rx('update job')
CC      call smvxt2(100+10*mod(job,10)+1,s_lat,s_pot,n1,n2,n3,1d-5,smrho)
CC      call smvxt2(100+10*mod(job,10)+0,s_lat,s_pot,n1,n2,n3,1d-3,smpot)
C
CC      call dbgsmrho('mkpot pot after smvxt2',k1,k2,k3,s_lat,smpot)
CC      call dbgsmrho('mkpot rho after smvxt2',k1,k2,k3,s_lat,smrho)
C#endif

C ... Add integrals of compensating gaussians with vext into gpot0, hpot0
C     gpot0-gpotb contributes to potential matrix element ppi (gaugm)
C     gpot0 made here should be exactly cancelled by gpotb (locpt2)
C     The corresponding energy contribution sgpe should cancel s_atparms%sgpote
C     Note that cv is not modified
      allocate(gpote(s_pot%nlml),hpote(nbas),cgsum(ng))
      call vesgcm(0,s_site,s_spec,s_lat,nbas,qmom,ng,
     .  cv,cgsum,f,gpote,hpote,qsmc,zsum,[s1])
      call prfrce(50,nbas,f,' Forces including (compensating gaussians) * vext:')
C     call yprm('gpote',1,gpote,0,s_pot%nlml,s_pot%nlml,1)
      call daxpy(s_pot%nlml,1d0,gpote,1,gpot0,1)
      call daxpy(nbas,1d0,gpote,1,hpot0,1)

C ... Back transform of density and potential to real-space mesh
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,1,0,1)
      call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)

C --- Collect energy terms; make zvnuc for smooth problem ---
      rcvxt = 0d0
      sgpe = 0d0
      iv0 = 0
      do  ib = 1, nbas
        is = s_site(ib)%spec
        call corprm(s_spec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        z = 0  ! Omit interaction of nuclear charge with vext
        lmxl = s_spec(is)%lmxl
        if (lmxl < 0) cycle
        nlm = (lmxl+1)**2
        rcvxt = rcvxt + (qcorg-z)*y0*gpote(iv0+1) + cofh*hpote(ib)
        do  ilm = 1, nlm
          sgpe = sgpe + qmom(iv0+ilm)*gpote(iv0+ilm)
        enddo
        iv0 = iv0+nlm
      enddo

      rhvxt = ue0 + sgpe

      deallocate(cv)
      deallocate(cgsum,gpote,hpote)

      call tcx('smvxt')

      end
C#ifdefC DBGVXT
C      subroutine smvxt2(job,s_lat,s_pot,n1,n2,n3,fac,smpot)
CC- Debugging printout routine for external potential
CC job
CC 100 digit 1 smpot is in k space
CC  10 digit 1 input, 0 output
CC   1 digit 1 smpot is density, not potential
C      implicit none
CC ... Passed parameters
C      integer job,n1,n2,n3
C      double precision fac
C      complex(8), target :: smpot(n1,n2,n3)
CC ... For structures
C      include 'structures.h'
C      type(str_lat)::   s_lat
C      type(str_pot)::   s_pot
CC ... Dynamically allocated local arrays
CC     complex(8), allocatable :: smfun(:,:,:)
C      complex(8), allocatable :: cv(:),cv0(:),cvn(:)
CC ... Local parameters
C      integer k1,k2,k3,ig,j,ng,ifi
C      double precision beta,xx,xx0,g2,tpiba,pi,pi8
C      procedure(integer) :: fxst,fopn
C      character fn*6
C      integer procid,mpipid
C      integer, parameter :: master=0
C
CC     return  ! Do nothing unless debugging
C
C      k1 = n1; k2 = n2; k3 = n3
C      ng = s_lat%ng
C      procid = mpipid(1)
C      pi   = 4d0*datan(1d0)
C      pi8  = 8*pi
C      tpiba=2*pi/s_lat%alat
C
C      if (procid /= master) goto 1
C
CC ... Switch to reciprocal space; collect coefficients
CC     allocate(smfun(n1,n2,n3))
CC     call dcopy(n1*n2*n3*3,smpot,1,smfun,1)
C      if (mod(job/100,10) == 0) call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,-1)
C      allocate(cv(ng))
C      call gvgetf(ng,1,s_lat%kv,k1,k2,k3,smpot,cv)
C
CC     If we are dealing with input potential or density, save it
C      if (mod(job/10,10) == 1) then
C        fn = 'smpotn'
C        if (mod(job,10) /= 0) fn = 'smrhon'
C        ifi = fopn(fn)
C        rewind ifi
C        write(ifi,"('% rows',i6)") ng
C        do  ig  = 1, ng
C          write(ifi,'(i6,2f18.12)') ig, cv(ig)
C        enddo
C        call fclose(ifi)
C        print "(/' wrote ',a)", fn
C      endif
C
CC     Read or write smpot0 or smrho0
C      fn = 'smpot0'
C      if (mod(job,10) == 1) fn = 'smrho0'
C
CC ... if no smpot file yet, create one and exit
C      if (fxst(fn) /= 1) then
C        ifi = fopn(fn)
C        rewind ifi
C        write(ifi,"('% rows',i6)") ng
C        do  ig  = 1, ng
C          write(ifi,'(i6,2f18.12)') ig, cv(ig)
C        enddo
C        call fclose(ifi)
C        print "(/' wrote ',a)", fn
C        call rx0('done')
C      endif
C
CC ... Print out changes in smpot
C      allocate(cv0(ng))
C      ifi = fopn(fn)
C      rewind ifi
C      print "(/' reading ',a,'; print changes relative to it')", fn
C      read(ifi,"(6x,i6)") ig
C      call sanrg(.true.,ig,ng,ng,'smvxt2','ng')
C      do  ig  = 1, ng
C        read(ifi,*) j, cv0(ig)
C      enddo
C      call fclose(ifi)
C
CC ... set cv(2) by hand
C      xx = cv(2); xx0 = xx
CC      if (mod(job/10,10) == 1 .and. mod(job,10) == 0) then
CC        print *, 'cv(2)', cv(2)
CC        write(*,"(/' Set to : ')",advance='no')
CC        read(*,*) xx
CC        cv(2) = xx; cv(3) = xx
CC      endif
C
C      if (mod(job,10) /= 0) then
C        call info0(1,0,0,'    ig       kv       pi8/q2%9frho%10frho-rho0%6fIm rho')
C      else
C        call info0(1,0,0,'    ig       kv       pi8/q2%9fphi%10fphi-phi0%6fIm phi')
C      endif
C      do  ig = 1, ng
CC       if (abs(s_pot%smvextc(ig)) > 1d-8) then
C        if (abs(cv(ig)-cv0(ig)) > fac) then
C          g2 = tpiba*tpiba*(s_lat%gv(ig,1)**2+s_lat%gv(ig,2)**2+s_lat%gv(ig,3)**2)
C          call info8(1,0,0,'%,6i%3,4i  %;12F  %;12F  %;12F  %;12F',ig,s_lat%kv(ig,:),
C     .      pi8/q2, cv(ig), cv(ig)-cv0(ig), dimag(cv(ig)), 7, 8)
C        endif
C      enddo
C
C      if (mod(job/10,10) == 1 .and. mod(job,10) == 0 .and. abs(xx-xx0) > 1d-8) then
C        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cv,smpot)
C      endif
C
C      if (mod(job/10,10) == 0 .and. mod(job,10) == 1) then
CC        Uncomment the following to restrict to local field corrections
CC        beta = .1d0
CC        beta = 1
CC        allocate(cvn(ng))
CC        fn = 'smrhon'
CC        ifi = fopn(fn)
CC        rewind ifi
CC        print "(/' massage density ... reading ',a)", fn
CC        read(ifi,"(6x,i6)") ig
CC        call sanrg(.true.,ig,ng,ng,'smvxt2','ng')
CC        do  ig  = 1, ng
CC          read(ifi,*) j, cvn(ig)
CC        enddo
CC        call fclose(ifi)
CC
CC        do  ig = 1, ng
CC          if (abs(s_pot%smvextc(ig)) > 1d-8) then
CC            g2 = tpiba*tpiba*(s_lat%gv(ig,1)**2+s_lat%gv(ig,2)**2+s_lat%gv(ig,3)**2)
CC            call info8(1,0,0,'%,6i%3,4i  %;12F  %;12F  %;12F  %;12F',ig,s_lat%kv(ig,:),
CC     .        pi8/q2, cv(ig), cv(ig)-cv0(ig), dimag(cv(ig)), 7, 8)
CCC           write(*,'(i6,3i4,2f15.10,2x,2f15.10,2x)') ig,s_lat%kv(ig,:), cv(ig), cv(ig)-cv0(ig)
CC            cv(ig) = beta * cv(ig) + (1-beta) * cvn(ig)
CC          else
CC            cv(ig) = cvn(ig)
CC          endif
CC        enddo
CC        call gvputf(ng,1,s_lat%kv,k1,k2,k3,cv,smpot)
CC        deallocate(cvn)
C      endif
C      deallocate(cv0)
C
C   99 continue
C      if (mod(job/100,10) == 0) call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
C      deallocate(cv)
C
C    1 continue
C      call mpibc1(smpot,2*size(smpot),4,.false.,'smvxt2','smpot')
C
C      end
C
C      subroutine dbgsmrho(msg,n1,n2,n3,s_lat,smrho)
CC- Convert smrho to k space for viewing components via a debugger
C      implicit none
CC ... Passed parameters
C      character msg*(*)
C      integer n1,n2,n3,ng
C      integer k1,k2,k3,kkk
CC ... For structures
C      include 'structures.h'
C      type(str_lat)::   s_lat
C
C      double complex smrho(n1,n2,n3)
C      complex(8), allocatable :: wk(:),cn(:)
C
C      k1 = n1; k2 = n2; k3 = n3
C      kkk = k1*k2*k3
C      ng = s_lat%ng
C
C      call info0(10,1,0,' dbgsmrho : '//msg)
C
C      allocate(cn(ng),wk(kkk))
C      call dcopy(kkk*2,smrho,1,wk,1)
C      call fftz3(wk,n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,s_lat%kv,k1,k2,k3,wk,cn)
C      deallocate(cn,wk)
C
C      end
C#endif
