      subroutine getw0l(procid,llw,iw1,iw2,nq0i,dmlx,epinvq0i,wklm,wbz,
     .lmxax,nphon,ntphon,wtl,freq,lreal,epinv,w0,llmat)
C- Effective screened Coulomb interaction W-V at q=0, with possible phonon correction
C ----------------------------------------------------------------------
Ci Inputs
Ci   procid: if 0, print results to stdout
Ci   llw   : <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
Ci   iw1,iw2 : w0 is calculated for frequency points iw1:iw2
Ci         : For frequencies on the real axis:
Ci         : iw1 = 0 when there is time reversal symmetry
Ci         : iw1 = -iw2 otherwise
Ci         : For frequencies on the imaginary axis:
Ci         : iw1 = 1
Ci   nq0i  : Number of offset-gamma points
Ci   dmlx  : coefficients to Ylm expansion of <ehat|epinv|ehat> where :
Ci   epinvq0i : <q0i/|q0i|| epinv(:,:,iq0j)|q0i/|q0i|> = \sum_{ij}^{} {q_i^m} epinv_{ij}^{m'} q_j^m /({q^m})^2
Ci         :    epinv is "invariant tensor" : symmetrized q_i q_j, q = offset point (diele_invariant)
Ci   wklm  : coefficient w_L for Auxiliary function in Eqns. 28,33 of JPSP 83, 094711 (2014) (getwklm)
Ci   wbz   : sampling weight of the central cell for BZ integration = 1/(n1*n2*n3)
Ci   lmxax : l-cutoff for K_lm
Ci   q0i   : offset Gamma k-points
Ci   ntphon: total number of phonons over all non-equivalent offset Gamma points
Ci             if not zero LST correction included
Ci  nphon(nq0i): number of phonons for each direction corresponding to offset Gamma points
Ci   wtl : wtl(1,:) = transverse optical frequency; wtl(2,:) = longitudinal optical frequency
Ci           for each mode for each inequivalent direction  (in Hartree)
Ci   freq  : freq(iw1:iw2) contain frequencies used in the Lyddane-Sachs-Teller factor
Ci         : They should correspond to w0(iw1:iw2).  Used only if wT0L0 is nonzero.
Ci         : Note: Without time reversal symmetry, iw1<0 on the real axis.  But freq is kept for w>=0 only.
Ci         : For freq(w<0) should be taken from
Ci   iwf1  : First frequency index for freq.
Ci   lreal : 1 if freq is real, 0 freq is imaginary, 2 if freq is complex (not implemented)
Co Output
Co   w0    : cell averaged effective screened coulomb interaction, singularity averaged
Co           Lyddane Sachs Teller correction included if ntphon /= 0 and wtl passed
Cr Remarks
Cr   W is singular as q->0, as 1/q^2.   The singularity is integrable.
Cr   This routine returns an effective averaged W at q=0, lhs of Eq. 33 of JPSP 83, 094711 (2014).
Cr
Cr   A correction from phonons of the Lyddane-Sachs-Teller type is included if wT0L0 is nonzero.
Cu Updates
Cu   17 Sep 14 A first cut at including a Lyddane-Sachs-Teller scaling of W(q=0). Adapted from getw0
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: procid,iw1,iw2,nq0i,lmxax,lreal,ntphon
      integer, intent(in) :: nphon(nq0i)
      real(8), intent(in) :: epinvq0i(nq0i,nq0i),dmlx(nq0i,9),wklm((lmxax+1)**2),wbz,epinv(3,3,nq0i),
     .  freq(max(iw1,0):iw2)
      complex(8),intent(inout):: llw(iw1:iw2,nq0i)
      complex(8),intent(out):: w0(iw1:iw2),llmat(3,3)
C ... Local parameters
      logical llds  ! T if incorporate Lyddane-Sachs-Teller correction
      integer lm1x,lm2x,iw,lm,lm1,lm2,nlklm
      real(8) :: fpi,pi,v0
      real(8) :: epinvq0i_m1(nq0i,nq0i)
      complex(8) :: omg,llwyl(9,iw1:iw2),llw_invr(iw1:iw2,1:nq0i),cgllw((lmxax+1)**2,(lmxax+1)**2)
C     real(8),allocatable:: cy(:),yl(:),yl2(:)
c      complex(8), allocatable :: w00(:)
      character*(80) strn
      integer :: iphon,ixx,iq0
!       real(8), intent(in) :: wtl(2,ntphon) ! issues with ntphon=0
      real(8), intent(in) :: wtl(2,*) ! issues with ntphon=0
      real(8) :: faclst, symops(9)
      real(8), allocatable :: cg(:,:,:)
#ifdef COMMONLL
      integer ll(51**2)
      common/llblock/ll
#else
      integer ll
      external ll
#endif

      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      v0 = wklm(1)*fpi**1.5d0/wbz
      nlklm = (lmxax+1)**2
      epinvq0i_m1 = epinvq0i
       llds= ntphon /= 0
C   renormalize <e|L|e> in llw for phonons by Lyddane-Sachs-Teller relation
      if (llds) then
         print *, 'Lyddane  Sachs Teller correction included'
      do iw=iw1,iw2
         if (iw < 0) call rx('check freq entering in getw0l')
              omg = dsign(freq(iabs(iw)),dble(iw)) ; if (lreal == 0)
     .        omg = dcmplx(0d0,dsign(freq(iabs(iw)),dble(iw)))
              if (iw < 0) omg = -omg
         ixx=0
         do iq0=1,nq0i
            faclst=1d0
C            print *, 'iq0=',iq0
            do iphon=1,nphon(iq0)
c               print *,'iphon=',iphon
            ixx=ixx+1
c            print *, 'omg2=',omg**2,'wt=',wtl(1,ixx),'wl=',wtl(2,ixx),
c     .' fac=',(wtl(2,ixx)**2-omg**2)/(wtl(1,ixx)**2-omg**2)
            faclst=faclst*(wtl(2,ixx)**2-omg**2)/(wtl(1,ixx)**2-omg**2)
            enddo
            llw(iw,iq0)=llw(iw,iq0)*faclst
            print *, 'LST factor=',faclst, 'for iq0=', iq0, 'iw=',iw,'frequency=',omg
         enddo
      enddo
      endif

      if (procid == 0) write(6,333) nq0i,iw1,iw2,lmxax,sum(abs(llw(iw1:iw2,1:nq0i)))
  333 format(' getw0l :',i2,' offset point(s)   Frequency mesh (',2i3,
     .  ')  Klm to lmx =',i2,'.  Sum check llw =',1pd13.5)
      call matinv(nq0i,epinvq0i_m1)

C ... llw_invr: representation of epsilon as a linear combination of invariate tensors
      llwyl = 0d0
      do iw = iw1, iw2
        llw_invr(iw,1:nq0i) = matmul (epinvq0i_m1, llw(iw,1:nq0i))
        lm = 1
        llwyl(lm,iw) = sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
        do lm = 5,9
          llwyl(lm,iw) = sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
        enddo
      enddo

!! omega=0 3x3 matrix, only when iw1<= 0 <= nw (only when including iw=0 for omega=0)
      if (iw1 <= 0) then
        iw = 0
        llmat = 0d0
        do  iq0 = 1, nq0i
          llmat(1:3,1:3) = llmat(1:3,1:3) + llw_invr(iw,iq0)*epinv(1:3,1:3,iq0)
        enddo
      endif

      allocate(cg(nlklm, nlklm, (2*lmxax+1)**2))
      cg = 0
      symops = [1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0]

      call rotcg(lmxax,symops,1,cg) !no rotation

C --- Loop over frequencies ---
      do  iw  = iw1, iw2
        lm2x = 0
        cgllw = 0d0
        do  lm2 = 1, nlklm
          if (mod(ll(lm2),2)==1) cycle  !only even l
          lm2x = lm2x+1
          lm1x = 0
          do  lm1 = 1, nlklm
            if (mod(ll(lm1),2)==1) cycle !only even l
            lm1x = lm1x+1
            cgllw(lm1x,lm2x) = cg(lm2,1,lm1)*llwyl(1,iw) + sum(cg(lm2,5:9,lm1)*llwyl(5:9,iw))
          enddo
        enddo

C   ... Make Klm (defined in PRB 81, 125102, Eq 43)
C       Inversion of 1 = \sum_L1 llwyl(L1)Y_L1 * \sum_L2 K_L2 Y_L2
C       Both sides are multiplied by Y_L and integrated over the sphere.
C       sqrt(fpi) comes from \int d\Omega Y_0 (right hand side of inversion eq).
        call matcinv(lm2x,cgllw(1:lm2x,1:lm2x))
        cgllw(1:lm2x,1:lm2x) = sqrt(fpi)*cgllw(1:lm2x,1:lm2x)

C   ... Spherical average of Klm
        lm1x = 0
        w0(iw) = 0
        do lm1 = 1, nlklm
          if (mod(ll(lm1),2)==1) cycle
          lm1x = lm1x+1
          w0(iw) = w0(iw) + wklm(lm1) * fpi * cgllw(lm1x,1)/wbz ! Klm=cgllw(:,1)
!         Takao Apr 2019 makes change ~/ecalj-master-from-takao-Apr2019/fpgw/gwsrc/mptauof.F
!         if(lm1==1) w0(iw)= w0(iw) - wklm(lm1)*fpi**1.5d0/wbz   ! fpi**1.5 --> subtract Coulomb
!         w0(iw)= wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).
        enddo

C   ... Lyddane-Sachs-Teller old version
C        w0(iw) = w00(iw)
C        if (llds) then
C          if (iw < 0) call rx('check freq entering in getw0l')
C          omg = dsign(freq(iabs(iw)),dble(iw)) ; if (lreal == 0) omg = dcmplx(0d0,dsign(freq(iabs(iw)),dble(iw)))
C          if (iw < 0) omg = -omg
C          w0(iw) = w00(iw) * (wT0L0(1)**2 - omg**2) / (wT0L0(2)**2 - omg**2)
CC          w0(iw) = w00(iw) * (wT0L0(1)**2 - 0*omg**2) / (wT0L0(2)**2 - 0*omg**2)
C        endif

C   ... Subtract Coulomb
        w0(iw) = w0(iw) - v0
c        w00(iw) = w00(iw) - v0
C       w00(iw) = wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).

      enddo ! frequency loop

      if (procid == 0) then
        strn = 'Re' ; if (lreal == 0) strn = 'Im'
        if (llds) then
          write(6,"(3x,a,' freq   iw',14x,'(W-v)0',18x,'(W-v)0 Modified by L-T-S')") trim(strn)
        else
          write(6,"(3x,a,' freq   iw',14x,'(W-v)0')") trim(strn)
        endif
        do  iw = iw1, iw2
          omg = dsign(freq(iabs(iw)),dble(iw))
c          if (llds) then
c            write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8)') dble(omg),iw,w00(iw),w0(iw)
c          else
            write(6,'(f10.4,i5,1p,2d16.8,2x,2d16.8)') dble(omg),iw,w0(iw)
c          endif
        enddo
        if (procid == 0) then
c          if (iw1 <= 0) write(6,334) v0, (dble(w00(0))+v0)/v0, (dble(w0(0))+v0)/v0, sum(abs(w0))
c  334     format(' Bare v0 =', 1pd13.5,'  W0/v0 =',1pd13.5,'  W0/v0(LDS) =',1pd13.5,'  Sum check W0 =',1pd13.5)
          if (iw1 <= 0) write(6,334) v0, v0*wbz, (dble(w0(0))+v0)/v0, sum(abs(w0))
  334     format(' Bare v0 =', 1pd13.5,' v0*vol0 =', 0pf11.5,'  W0/v0 =',1pd13.5,'  Sum check W0 =',1pd13.5)
          if (iw1 > 0) write(6,"(' Bare v0 =', 1pd13.5,'  Sum check W0 =',1pd13.5)")
     .    v0, sum(abs(w0))
        endif
      endif
c      deallocate(w00)

      deallocate(cg)
      end

C      subroutine getw0(llw,ii,ie,nq0i,dmlx,epinvq0i,wklm,wbz,lmxax,q0i,epinv, w0)
CC- Obtain effective screened Coulomb interaction W-V at q=0
CCi Inputs
CCi   llw   : <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
CCi   ii    : w0 is calculated for points ii:ie
CCi   ie    : w0 is calculated for points ii:ie
CCi   nq0i  : Number of offset-gamma points
CCi   dmlx  : coefficients Ylm expansion of <ehat|epinv|ehat> where
CCi   epinv : is "invariant tensor" : symmetrized q_i q_j, q = offset point
CCi   epinvq0i : \sum_{ij}^{} {q_i^m} epsinv_{ij}^{m'} q_j^m /({q^m})^2
CCi   wklm  : related to K_lm defined around Eq.35 in Computer Physics Comm. 176,1-13 (2007)
CCi         : or maybe rather coefficient w_L in Eqns. 28,33 of JPSP 83, 094711 (2014)
CCi   wbz   : sampling weight for BZ integration
CCi   lmxax : L-cutoff for K_lm
CCi   q0i   : offset Gamma k-points
C!! Output
C!!   w0(ii:ie), others are inputs.
Cu Updates
Cu    Superseded by getw0l
CC ----------------------------------------------------------------------
C      implicit none
C      integer:: ii,ie,nq0i,lm,lm1,lm2,nlklm,lmxax,iw,ll,lm1x,lm2x
C      complex(8):: llw(ii:ie,nq0i),w0(ii:ie)
C      complex(8):: llwyl(9,ii:ie),llw_invr(ii:ie,1:nq0i)
C      real(8):: epinvq0i(nq0i,nq0i),epinvq0i_m1(nq0i,nq0i),dmlx(nq0i,9)
C      real(8)::  cg((lmxax+1)**2,(lmxax+1)**2,(2*lmxax+1)**2)
C      complex(8)::   cgllw((lmxax+1)**2,(lmxax+1)**2),sss,sss2
C      real(8):: wklm((lmxax+1)**2),ylm,wbz,rrr2(3)
C      real(8):: fpi,pi,epinv(3,3,nq0i),q0i(3,nq0i),  r2s,emat(3,3),rrr(3),qnorm2(nq0i)
C      integer:: nlmm,lx,iq0i,iq0x
C      real(8),allocatable:: cy(:),yl(:),yl2(:)
C
C      external ll
C
C
C      pi  = 4d0*datan(1d0)
C      fpi = 4d0*pi
C      nlklm=(lmxax+1)**2
C      epinvq0i_m1=epinvq0i
C      write(6,"(' getw0 start:',4i3,d13.6)") nq0i,lmxax,ii,ie,sum(abs(llw(ii:ie,1:nq0i)))
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc      nlm=(lmxax+1)**2
Cc      do lm=1,nlm
Cc        if(abs(wklm(lm))>1d-10) write(6,*)lm,ll(lm),wklm(lm)
Cc      enddo
Cc      stop 'test end vvvvvvvvvvvvvvvvvvvvvv'
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C      call matinv(nq0i,epinvq0i_m1)
C      llwyl=0d0
C      do iw= ii,ie
C        llw_invr(iw,1:nq0i) = matmul (epinvq0i_m1, llw(iw,1:nq0i))
C        ! llw_invr: representation of epsilon as a linear combination of invariate tensors
C        lm=1
C        llwyl(lm,iw) = sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
C        do lm=5,9
C          llwyl(lm,iw)= sum(dmlx(1:nq0i,lm)*llw_invr(iw,1:nq0i))
C        enddo
Cc   write(*,"('qqq iw llw',i3,3(2d12.4,x),3x,20(2d12.4,2x))")iw,llw(iw,1:nq0i),llwyl(1,iw),llwyl(5:9,iw)
C      enddo
C      call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg) !no rotation
C
Cc$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc$$$      lx=2
Cc$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
Cc$$$      call sylmnc(cy,lx)
Cc$$$      do iw= ii,ie
Cc$$$        do iq0x=1,nq0i
Cc$$$        sss=0d0
Cc$$$        do iq0i=1,nq0i
Cc$$$            sss= sss  + llw_invr(iw,iq0i)*sum(q0i(:,iq0x)*matmul(epinv(:,:,iq0i),q0i(:,iq0x)))
Cc$$$     &                 /sum(q0i(:,iq0x)**2)
Cc$$$        enddo
Cc$$$          write(*,"(' ttt: epinv expansion=',2i4,2f10.5,2x,2d13.5)") iq0x,iw,sss,sss-llw(iw,iq0x)
Cc$$$        enddo
Cc$$$        write(*,*)
Cc$$$      enddo
Cc$$$      stop '---- ttt test1: reproduce llw at q0i -----------------------'
C
Cc$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc$$$!! === test for one r vector as for <ehat| epinv|ehat> = \sum_lm dmlx(iq0i,lm) *Y_lm(ehat) ===
Cc$$$!! ===== generate YL for a test vector rrr (rrr is ehat above).=====
Cc$$$      lx=2
Cc$$$      allocate(cy((lx+1)**2),yl((lx+1)**2),yl2((lx+1)**2))
Cc$$$      call sylmnc(cy,lx)
Cc$$$      do i=1,3
Cc$$$         if(i==1) rrr =(/1d0,1d0,0d0/)
Cc$$$         if(i==2) rrr =(/1d0,-1d0,0d0/)
Cc$$$         if(i==3) rrr =(/0d0,0d0,1d0/)
Cc$$$         rrr = rrr/sqrt(sum(rrr**2))
Cc$$$c     write(*,"(' testttt: r=',3f10.5)") rrr
Cc$$$         call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
Cc$$$         do iw= ii,ie
Cc$$$            sss=0d0
Cc$$$            do iq0i=1,nq0i
Cc$$$               sss = sss + llw_invr(iw,iq0i)*sum(rrr*matmul(epinv(:,:,iq0i),rrr))
Cc$$$            enddo
Cc$$$c            if(abs(llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss)>1d-12) then
Cc$$$            write(*,"(' ttt: epinv expansion=',i3,2f10.5,2x,4d13.5)") iw,sss,
Cc$$$     &      llwyl(1,iw)*cy(1)*yl(1)+sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9))-sss, llw(iw,i)-sss
Cc$$$c            endif
Cc$$$         enddo
Cc$$$      enddo
Cc$$$      stop '---- ttt testxxx: reproduce llw at q0i -----------------------'
C
Cc$$$      lx=lmxax
Cc$$$      if(allocated(cy)) deallocate(cy)
Cc$$$      if(allocated(yl)) deallocate(yl)
Cc$$$      allocate(cy((lx+1)**2),yl((lx+1)**2))
Cc$$$      call sylmnc(cy,lx)
Cc$$$      rrr=(/-0.36d0,0.20d0,0.4d0/)
Cc$$$      rrr=rrr/sqrt(sum(rrr**2))
Cc$$$      call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
Cc$$$      write(*,"(  ' test input: q=',3f10.5)") rrr
C
C!! ---
C      w0=0d0
C      do iw =ii,ie
C        lm2x=0
C        cgllw=0d0
C        do lm2=1,nlklm
C          if(mod(ll(lm2),2)==1) cycle  !only even l
C          lm2x=lm2x+1
C          lm1x=0
C          do lm1=1,nlklm
C            if(mod(ll(lm1),2)==1) cycle !only even l
C            lm1x=lm1x+1
C            cgllw(lm1x,lm2x) = cg(lm2,1,lm1)*llwyl(1,iw) + sum(cg(lm2,5:9,lm1)*llwyl(5:9,iw))
C          enddo
C        enddo
C
C!! === inversion of 1 = \sum_L1 llwyl(L1)Y_L1 * \sum_L2 K_L2 Y_L2===
C!! Both sides are multipled by Y_L and integrated over the sphere.
C!! warn: Klm should be not confused with wklm
C!! Klm is defined in Christoph's paper PRB 81, 125102-8 Eq 44.
C        call matcinv(lm2x,cgllw(1:lm2x,1:lm2x))
C        cgllw(1:lm2x,1:lm2x)= sqrt(fpi)*cgllw(1:lm2x,1:lm2x)
C        !sqrt(fpi) comes from \int d\Omega Y_0 (right hand side of inversion eq).
C
C!! === spherical integral of Klm ===
C        lm1x=0
C        do lm1=1, nlklm
C          if(mod(ll(lm1),2)==1) cycle
C          lm1x=lm1x+1
C          w0(iw)= w0(iw)+ wklm(lm1)* fpi* cgllw(lm1x,1)/wbz      ! Klm=cgllw(:,1)
C          if(lm1==1) w0(iw)= w0(iw) - wklm(lm1)*fpi**1.5d0/wbz   ! fpi**1.5 --> subtract Coulomb
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc          w0(iw)= wklm(1)*fpi**1.5*(1d0/llw(iw,1)-1d0)/wbz !test case of llw at q0i(:,iq0i=1).
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C        enddo
C
Cc$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cc$$$        lx=lmxax
Cc$$$        print *,' lx=',lx
Cc$$$        if(allocated(cy)) deallocate(cy)
Cc$$$        if(allocated(yl)) deallocate(yl)
Cc$$$        allocate(cy((lx+1)**2),yl((lx+1)**2))
Cc$$$        call sylmnc(cy,lx)
Cc$$$        rrr=(/1d0,1d0,0d0/)
Cc$$$c        rrr=(/1d0,-1d0,0d0/)
Cc$$$c        rrr=(/0d0,0d0,1d0/)
Cc$$$        rrr=rrr/sqrt(sum(rrr**2))
Cc$$$        call sylm(rrr,yl,lx,r2s) !spherical factor Y( q+G )
Cc$$$c        write(*,"(  ' qqq test input: q=',3f10.5)") rrr
Cc$$$        sss=0d0
Cc$$$        lm2x=0
Cc$$$        do lm2=1,nlklm
Cc$$$           if(mod(ll(lm2),2)==1) cycle !only even l
Cc$$$           lm2x=lm2x+1
Cc$$$           sss= sss + cgllw(lm2x,1) *cy(lm2)*yl(lm2)
Cc$$$        enddo
Cc$$$        write(*,"(' qqq ep test:',i3,2f10.5,2x,2f10.5,2x,2f10.5)") iw,sss,
Cc$$$     &       1d0/llw(iw,1)
Cc$$$c     &      1d0/(llwyl(1,iw)*cy(1)*yl(1) + sum(llwyl(5:9,iw)*cy(5:9)*yl(5:9)))
Cc$$$      stop ' --- test3 ttt vvvvvvvvvvvvvvvvv test end vvvvvvvvvvvvvvvvv ---'
C      enddo
Cc$$$      do iw=ii,ie
Cc$$$        write(6,"('w0 & w0(1:nq0i)=',2d11.3,10(1x,2d11.3))")
Cc$$$     &  w0(iw), (wklm(1)*fpi**1.5*(1d0/llw(iw,iq0i)-1d0)/wbz,iq0i=1,nq0i)
Cc$$$      enddo
C      end

      subroutine eibzgen(nqibz,symgg,ngrp,qibze,iqxini,iqxend,qbz,nqbz,timereversal,ginv,iprintx,
     o nwgt,igx,igxt,eibzsym,timerout)
!! === Obtain info for eibz symmetrization. See PRB81,125102-9 ===
!! For GaAs 4x4x4 (with timereversal), we have IBZxBZ=10x640 is reduced to be IBZxEBZ=286.
      implicit none

      integer:: nqibz,ngrp,iqxini,iqxend,iq,ig,neibzx,ibz,nwgtsum,ntimer
      integer:: neibz(iqxini:iqxend),nwgt(nqbz,iqxini:iqxend),ik,nqbz,i,it
      logical:: timereversal
      real(8):: symgg(3,3,ngrp),qibze(3,iqxini:iqxend),qbz(3,nqbz),
C     &    qeibz(3,nqbz,iqxini:iqxend),
     &  q(3),qx(3),qlat(3,3),qdiff(3),ginv(3,3),qxx(3),timer,qxi(3) ,
     &  sss(3,3),sumcheck1 !qout(3,nqbz,iqxini:iqxend,nqbz)
      integer:: igx(ngrp*2,nqbz,iqxini:iqxend),igxt(ngrp*2,nqbz,iqxini:iqxend),imm,im
      integer:: eibzsym(ngrp,-1:1,iqxini:iqxend),ibzx,eallow(ngrp,-1:1),immxx
      logical:: iprintx,timerout

      integer,allocatable:: immx(:),mqbz(:)
      integer,allocatable:: nxx(:,:,:),igx_(:,:),igxt_(:,:)
C     integer,allocatable:: ibznxx(:,:)

      real(8)::    epd=1d-5
      integer:: nkey(3),isig,kkk,kkk3(3),ik1,ik2,ik3,ix,ifi,pwmode,ifiese
      integer,allocatable:: ieord(:)
      integer,allocatable:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      real(8),allocatable:: qbzl(:,:)
      logical :: inotime,nexist
!!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cnnn
      real(8):: alatx, plat(3,3)
      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx)
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3)
      ifi = 6661
      open (ifi, file='LATTC')
      read(ifi,*) alatx
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      close(ifi)
      call minv33tp(plat,qlat)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(iprintx) write(6,"('eibzgen: TimeReversalSwitch ngrp= ',l,i3)") timereversal,ngrp
      ntimer=1
      if(timereversal) ntimer=-1
      sss(1,:)=(/1d0,0d0,0d0/)
      sss(2,:)=(/0d0,1d0,0d0/)
      sss(3,:)=(/0d0,0d0,1d0/)
      sumcheck1= sum(abs(symgg(:,:,1)-sss))
Cstop2rx 2013.08.09 kino      if(sumcheck1>1d-8) stop 'eibzgen: symgg(:,:,1) is not E'
      if(sumcheck1>1d-8) call rx( 'eibzgen: symgg(:,:,1) is not E')

!!=== Get key and nkey for each ix. (similar with those in readqg,readeigen) ===
      allocate(qbzl(3,nqbz),key(3,0:nqbz),ieord(nqbz))
      key(:,0)=0 !dummy
      do ibz=1,nqbz
        call rangedq(matmul(ginv,qbz(:,ibz)), qbzl(1,ibz))
      enddo
      do ix =1,3
         call sortea(qbzl(ix,:),ieord,nqbz,isig)
         ik=0
         do i=1,nqbz
            kkk= ( qbzl(ix,ieord(i))+0.5*epd )/epd  !kkk is digitized by 1/epd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
c               write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
Cstop2rx 2013.08.09 kino               stop 'iqindx: bug not sorted well'
               call rx( 'iqindx: bug not sorted well')
            endif
         enddo
         nkey(ix)=ik
      enddo
      deallocate(ieord)
!!  key is reallocated. inverse mattping, iqkkk
      allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
      kk1(:) = key(1,1:nkey(1))
      kk2(:) = key(2,1:nkey(2))
      kk3(:) = key(3,1:nkey(3))
      deallocate(key)
      allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
c      print *,' initqqq nqtt=',nqtt
      do i=1,nqbz
         kkk3= (qbzl(:,i)+0.5*epd)/epd !kkk is digitized by 1/epd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
c         print *, ' ik1,ik2,ik3,iqkkk=',ik1,ik2,ik3,i
      enddo
      deallocate(qbzl)
!! -----------------------------------------------------------
ccccccccccccccccccc
C      ifiese=9039
C      inquire(unit=ifiese,opened=nexist)
C      if(nexist) stop 'eibzgen: I assumed ifio=9039 is not used(too sloppy) pwmode'
      pwmode = 0
      INQUIRE(FILE='PWMODE',EXIST=nexist)
      if(nexist) then
        ifiese=9039
        open(ifiese,file='PWMODE')
        read(ifiese,*,err=999,end=999) pwmode
        close(ifiese)
      endif
cccccccccccccccccc
      if(pwmode>0.and.pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)

!! ===  main iq loop ===
      eibzsym = 0
      do iq = iqxini,iqxend
c        write(6,"('iq =',i7,'  xxxxxxxxxxxxxxxxxxxxxxx')") iq
        q = qibze(:,iq) !q means k in eq.50 of PRB81,125102-9
!! Allowed operation to keep q,   eallow(ig,it) \in EIBZ(q)
        eallow=0d0
        do it=1,ntimer,-2
          timer=dble(it)
        do ig=1,ngrp
          qdiff= timer*matmul(symgg(:,:,ig),q) - q
          call rangedq(matmul(ginv,qdiff), qxx)
c        if(iprintx) print *,'ig qdiff',ig,qdiff,q,ginv,qxx
cccccccccccccccccccccccccccccccccccccccccccccc
          if(pwmode>0.and.pwmode<10) then
c     check whether q is on 1st BZ boundary or not
             ppin=matmul(transpose(plat),q)
             call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
             if(ig/=1.and.nout>1) then
                do iout=1,nout
                   write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout),
     &                  sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &                  matmul(qlat(:,:),ppin+nlatout(:,iout))
                enddo
                print *
                cycle
             endif
          endif
ccccccccccccccccccccccccccccccccccccccccccccc
          if(sum(abs(qxx))<1d-8) then
            eallow(ig,it)=1
          endif
        enddo
        enddo
        eibzsym(:,:,iq)=eallow ! === eibzsym is eallow -----

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        print *,'nqbz=',nqbz,ngrp,ntimer, 8*nqbz*ngrp/1000000000.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! === Generate all nxx (integered q with epditit) generated from ibz. ===
        allocate(nxx(3,nqbz,ngrp*2),immx(nqbz))
        allocate(igx_(ngrp*2,nqbz),igxt_(ngrp*2,nqbz))
        immxx=0
        do ibz=1,nqbz
          imm=0
          do it=1,ntimer,-2 !ntimer is 1 or -1.
            timer = dble(it)
          do ig=1,ngrp
            if(eallow(ig,it)/=1) cycle !eallow=1 keeps q invariant. Allowed operation of EIBZ(q)
            qx = timer*matmul(symgg(:,:,ig),qbz(:,ibz))
            qxi= matmul(ginv,qx)
            call rangedq(qxi, qxx) ! qxx(i) takes from 0 to 1.
            do im=1,imm
              if(sum(abs(nxx(:,ibz,im) - (qxx(:)+.5*epd)/epd)) <2) then
                 goto 3012      !<2 for safe
              endif
            enddo
            imm=imm+1
            nxx(:,ibz,imm)= (qxx+.5*epd)/epd  !integer nxx(i) takes from 0 to 1/epd.
            igx_ (imm,ibz) = ig  !inefficient memory usage, I think.
            igxt_(imm,ibz) = it  !
c            write(6,"(' imm ibz it ig  nxx=',4i7,'  ',3i6)")imm,ibz,ig,it ,nxx(:,ibz,imm)
!! inversion table:
!!             ibz<---from qbzrid
!! we can skip such ibz.
 3012       continue
          enddo
          enddo
          immx(ibz)=imm
          if(imm>immxx) immxx=imm
        enddo

c        print *,'eibz aaaaaaa'
c        call cputid(0)

!! inequivalent points.
          allocate(mqbz(nqbz))
          inotime=.true.
          do ix=1,ntimer,-2 !try no time reversal for ix=1
             mqbz=0
             nwgt(:,iq)=0
             igx(:,:,iq) =999991
             igxt(:,:,iq)=999992
             do ibz=1,nqbz
                if(mqbz(ibz)==1) cycle
                do im=1,immx(ibz)
                   if(inotime) then
                      if(igxt_(im,ibz) ==-1) cycle
                   endif
                   kkk3 = nxx(:,ibz,im)
c                   print *,' ibz,im kkk3=',ibz,im,kkk3
                   call tabkk(kkk3(1), kk1,nkey(1), ik1)
                   call tabkk(kkk3(2), kk2,nkey(2), ik2)
                   call tabkk(kkk3(3), kk3,nkey(3), ik3)
                   ibzx=iqkkk(ik1,ik2,ik3)
                   if(mqbz(ibzx)==1) then
                      cycle
                   else
                      mqbz(ibzx)=1
                      nwgt(ibz,iq)= nwgt(ibz,iq) + 1
                      igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)
                      igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)
                   endif
                enddo
             enddo
             if(sum(nwgt(:,iq))==nqbz) then
                if(ix==1) eibzsym(:,-1,iq)=0
                goto 1202
             endif
             inotime=.false.
             print *,' mptauo: we need time reversal'
          enddo
          call rx( 'eibzgen: bug. nwgt sum is not nqbz')
 1202     continue

C cccccccccccccccccccccccccccccccccc
C         do ibz=1,200 !nqbz
C           if(nwgt(ibz,iq)/=0) then
C             write(6,"('yyy0: ',i8,2x,25(i3,i2))") ibz,(igx(i,ibz,iq),igxt(i,ibz,iq),i=1,nwgt(ibz,iq))
C           endif
C         enddo
C cccccccccccccccccccccccccccccccccc

c          print *,'pppppp nwgt(ibz,iq),nqbz=',iq,sum(nwgt(:,iq)),nqbz
c$$$!! === Search other ibz generated from ibz. ===
c$$$!!   ibznxx(ibz,imm) means "ibz rotated by igx(imm,ibz,iq),igxt(imm,ibz,iq)"
c$$$        allocate(mqbz(nqbz))
c$$$        do ibzx=1,nqbz
c$$$          do ibz=1,ibzx !search equivalent (ibzx,1) with first appeard (ibz,im) .
c$$$          do im=1,immx(ibz)
c$$$             if(sum(abs(nxx(:,ibzx,1) - nxx(:,ibz,im))) <2) then
c$$$               nwgt(ibz,iq) = nwgt(ibz,iq) +  1
c$$$               igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)
c$$$               igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)
c$$$               goto 3312
c$$$             endif
c$$$          enddo
c$$$          enddo
c$$$ 3312     continue
c$$$        enddo
        deallocate(igx_,igxt_)

c        print *,'eibz bbbbbb'
c        call cputid(0)
!
        neibzx=0
        do ibz=1,nqbz
           if(nwgt(ibz,iq) /=0) then
              neibzx=neibzx+1
c              qeibz(:,neibzx,iq)=qbz(:,ibz)
c              write(6,"('   ibz qeibz =',i7,3f10.5, ' nwgt=',i7)")ibz,qeibz(:,neibzx,iq),nwgt(ibz,iq)
           endif
        enddo
        neibz(iq) = neibzx !number of eibz for iq
        if( iprintx ) then
          write(6,"('iq=',i8,' # of EIBZ: Full=',i8,
     &      ' Used(TimeR 1 or -1)=',i3,'=',i2,'+',i2,' neibz= ',i7)")iq,sum(abs(eallow(:,:))),
     &      sum(eibzsym(:,:,iq)),sum(eibzsym(:,1,iq)),sum(eibzsym(:,-1,iq)),neibz(iq)
c          write(6,"('eibz: iq neibz nqbz= ',i3,3f11.5,3i7)") iq,q,neibz(iq),nqbz
        endif
        deallocate(nxx,immx,mqbz)
      enddo
!!
      nwgtsum = sum(nwgt(1:nqbz,iqxini:iqxend))
      timerout=timereversal
      if(minval(igxt)==1) timerout=.false.

      if(iprintx) then
        write(6,"(' nqbz,  sum(neibz(iq)), sum(ngwt)=sum(nqbz)= ',10i7)")
     &  nqbz,sum(neibz(iqxini:iqxend)), nwgtsum, nqbz*(iqxend-iqxini+1)
      endif
c      stop 'test stop -------- end of eibzgen ----------'

      return
  999 continue
      call rx('file PWMODE is empty')
      end subroutine eibzgen
!! -------------------------------------------------
      subroutine tabkk(kkin, kktable,n, nout)
      integer:: nout,n, kkin, kktable(n),i,mm,i1,i2
      i1=1
      i2=n
      if(kkin==kktable(1)) then
         nout=1
         return
      elseif(kkin==kktable(n)) then
         nout=n
         return
      endif
      do i=1,n
         mm=(i1+i2)/2
         if(kkin==kktable(mm)) then
            nout=mm
            return
         elseif(kkin>kktable(mm)) then
            i1=mm
         else
            i2=mm
         endif
      enddo
c$$$      do i=1,n
c$$$         if(kkin==kktable(i)) then
c$$$            nout=i
c$$$            return
c$$$         endif
c$$$      enddo
      print *, i1,i2,kkin
      print *, kktable(i1),kktable(i2)
Cstop2rx 2013.08.09 kino      stop 'takk: error'
      call rx( 'takk: error')
      end subroutine tabkk



      subroutine mptauof(symops,ng,plat,nbas,bas,iclass,miat,tiat,invg,delta)
C- Mapping of each atom by point group operations
C Modded by okuda 1994 March.
C Simplified by kotani 1994 7/31.
Ci  Input
Ci     symops(1,ng),ng,plat,nbas,bas(3,nbas)
Ci     iclass(nbas); denote class for each atom
Co  Output
Co    miat(ibas,ig) = is atom that ibas is rotated into, though shifted by tiat.
Co    tiat(ibas,ig) = translation to add bas(:,miat(ibas,ig)) to coincide with bas', Eq (1).
Co    delta = shifting vector for non-symmorphic group, called ag in lm package
Cr  Remarks
Co    space group operation maps: r' = matmul (am, r) + delta
Cr
Cr    (1) Specifically the ibas-th atom (position at bas(k,ibas) ) is mapped to
Cr        bas' = bas(:,miat(ibas,ig) )+ tiat(:,ibas,ig)
Cr             = matmul(symop(ig),bas(:,ibas)) + delta(:,ig)
Cr
Cr    (2) tiat = shift to add to bas(:,miat(ibas,ig)) to coincide with bas', Eq (1).
Cr        tiat must be some integer multiple of plat
C--------------------------------------------------------------------
      use mpi_mod, rank => mpi_rank

      implicit none
      integer ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng),ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &                 tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &                 ep,dd1,dd2,dd3,t1,t2,t3
      integer  iprint
      external iprint
      integer ires(3, nbas, ng)
      integer(4):: ib1,ib2
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      if(rank == 0) write(6,*)' MPTAUO2: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ) <= ep.and.
     &        abs( symops(2,ig)-symops(4,igd) ) <= ep.and.
     &        abs( symops(3,ig)-symops(7,igd) ) <= ep.and.
     &        abs( symops(4,ig)-symops(2,igd) ) <= ep.and.
     &        abs( symops(5,ig)-symops(5,igd) ) <= ep.and.
     &        abs( symops(6,ig)-symops(8,igd) ) <= ep.and.
     &        abs( symops(7,ig)-symops(3,igd) ) <= ep.and.
     &        abs( symops(8,ig)-symops(6,igd) ) <= ep.and.
     &        abs( symops(9,ig)-symops(9,igd) ) <= ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
        if(iprint() >= 40) then
          if (rank == 0) then
            print *,' '
            print *,' '
            print *,' **** group ops no. ig (igd)= ', ig, invg(ig)
            write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
            write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
            write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
 1731       format (' ',3f9.4)
          endif
        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c       Look for translation vector ag (it is not supplied but calculated here and called tran
        do 20 ib1=1,nbas
        do 20 ib2=1,nbas
          tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))

c         (b1,b2,b3) = coordinates of ibas after rotation + shft by trial tran
          do 30 ibas=1,nbas
            b1=am(1,1)*bas(1,ibas)
     &        +am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas)
     &        +tran(1)
            b2=am(2,1)*bas(1,ibas)
     &        +am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas)
     &        +tran(2)
            b3=am(3,1)*bas(1,ibas)
     &        +am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas)
     &        +tran(3)
c
c           Find for atom mi for which (b1,b2,b3) - bas(:,mi) = lattice vector
            do 40 mi=1,nbas
              if( iclass(mi) /= iclass(ibas) ) go to 40

              do 50 i1=-3,3
              do 50 i2=-3,3
              do 50 i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1) <= ep.and.abs(t2) <= ep.and.
     &             abs(t3) <= ep) go to 60
   50         continue
   40       continue
C           Search failed. Try next (tr).
            goto 20

   60       continue
            miat(ibas,ig)  = mi
            tiat(1,ibas,ig)= dd1
            tiat(2,ibas,ig)= dd2
            tiat(3,ibas,ig)= dd3
            ires(1,ibas,ig)= i1
            ires(2,ibas,ig)= i2
            ires(3,ibas,ig)= i3
c
   30     continue
c When the do-30 loop has been completed, we get out of do-20 loop
          goto 21
   20   continue
        call rx( 'mptauo2: Can not find miat and tiat')
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

c- have gotten the translation-> check write --------------------
        if(iprint() >= 40) then

          if (rank == 0) then
            write(6,4658)tran
 4658       format('  Obtained translation operation=',3d12.4)
            do 123  ibas=1,nbas
              write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &      tiat(2,ibas,ig), tiat(3,ibas,ig),
     &      ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150         format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123       continue
          endif

C          print 554, 'check formula'//
C     .      ' bas(miat(ibas,ig)) = matmul(symop(ig),bas(ibas)) + delta(ig) - tiat(ibas,ig)'
C  554     format(a)
C          do  ibas=1,nbas
C            tran = matmul(am,bas(:,ibas)) + delta(:,ig) - tiat(:,ibas,ig)
C            tran = tran - bas(:,miat(ibas,ig))  ! This should be zero
C            print 555, ibas, tran+bas(:,miat(ibas,ig)), tran
C  555       format(i4,3f12.6,2x,3f12.6)
C            if (sqrt(tran(1)**2+tran(2)**2+tran(3)**2) .gt. 1d-6) stop 'oops'
C          enddo

        endif
   10 continue

      end
