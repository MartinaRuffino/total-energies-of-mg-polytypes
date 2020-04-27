      subroutine mksecp(iSigma,npm,nstate,ntq,ntqx,nt0p,nt0m,nt_max,
     .  nctot,ngb,ekc,ef,omega,GaussSmear,esmr,nw,wtkp,nwxi,nwx,nw_i,
     .  freq_r,zw3,zmel,zsec)
C- Adds pole contribution to self energy for QSGW Vxc
C ----------------------------------------------------------------------
Ci Inputs
Ci   iSigma: mode determining form of static approximation to SE
Ci         : 0 SE_nn'(ef)+imag integral:delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ci         : 1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
Ci         : 2 SE_nn'((e_n+e_n')/2)
Ci         : 3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
Ci         : 4 SE_nn'(ef)
Ci         : 5 delta_nn' SE_nn(e_n)
Ci   npm   :1 Use time reversal symmetry
Ci         :2 No time reversal symmetry
Ci   nstate:number of states contributing to sigma
Ci   ntq   :dimension of 1-particle hamiltonian; number of eigenstates t
Ci         :dimensions zsec
Ci   ntqx  :Generate SE for ntqx states
Ci   nt0p  :upper bound to it index when omg(itp)<ef
Ci   nt0m  :lower bound to it index when omg(itp)>ef
Ci   nt_max:upper bound to it index when omg(itp)>ef
Ci   nctot :number of core states
Ci   ngb   :Number of product basis functions
Ci   ekc   :ekc(1:nstate) QP levels
Ci   ef    :Fermi level
Ci   omega :frequencies for each QP level
Ci   GaussSmear
Ci         : T, Broaden poles in G with Gaussian of width esmr
Ci         : F, Broaden poles in G with rectangle of width esmr
Ci   esmr  :Broadening of G in integration; see wintzg
Ci   nw    :Number of imaginary mesh points
Ci   nw_i  :lower bound of freq_r
Ci   nwxi  :dimensions zw3 - zw3(:,:,nwxi:nwx)
Ci   nwx   :dimensions zw3 - zw3(:,:,nwxi:nwx)
Ci   wtkp  :k-point weight
Ci   freq_r:Frequencies on the real axis
Ci   zmel  :matrix elments of the product basis
Ci   zw3   :
Co Outputs
Co   zsec  :Pole contribution to S_ij= <i|Re[S](e_i)|j> for a given kp
Cs Command-line switches
Cl Local variables
Cl   iir   : 1 if current omega is > ef
Cl         :-1 if current omega is < ef
Cr Remarks
Cr  This contribution to zsec_nn' is not hermitian and must be made hermitian later,
Cr  as described in Remarks, routine mkseci
Cr  See mksecp for a vectorized version of this routine.
Cr
Cr  This routine has similarities to matzwz3.
Cr  We need to calculate a quantity like this:
Cr    zwz(t,t',t'',w) = S[I,J] zm+(I,t,t') zw3(I,J,w'[t,t']) zm(J,t,t'')
Cr  where
Cr    I,J are mixed basis indices
Cr    w is an index to an energy point which ranges from 1 to 3 (for interpolation)
Cr    w' is a related index pointing to the appropriate point in zw3
Cr    t,t',t'' are band indices.
Cr    zm is a subset of zmel, with some elements t omitted; see below.
Cr
Cr  We can structure the calculation this way:
Cr    zwz(t,t',t'',w )= S[J] CC(J,t,t',w) zm(J,t,t'')        (A)
Cr    where
Cr      CC(J,t,t',w) = S[I] zm+(I,t,t') zw3(I,J,w'[t,t'])   (B)
Cr    In practice this is
Cr      cc(t,t',w)_J,t,w
Cr
Cr  We cannot follow the logic of matzwz3 for the following reason:
Cr  There the (t,t') loop order was reversed:
Cr        do  t = 1, nstate
Cr          do t' = 1, ntqx
Cr            do t'' = 1, ntqx
Cr              zwzi(t)_t',t'' = sum_J CC(J,t,t') zm(J,t,t'')
Cr            enddo
Cr          enddo
Cr        enddo
Cr
Cr  The equivalent is not possible here, because
Cr  (1) only a portion of the states in loop [t = 1, nstate]
Cr  contribute to the result.
Cr  Which elements t are omitted depends on t'.
Cr
Cr  (2) Potential zw3 depends on energy through both t and t'.
Cr  The best we can do are make inner products through zgemv.
Cr
Cu Updates
Cu   22 May 13 Adapted from Takao's sxcf_fal2.sc.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical GaussSmear
      integer iSigma,npm,nstate,ntq,ntqx,nt0m,nt0p,nt_max,nctot,ngb,nw,nwxi,nwx,nw_i
      real(8) :: omega(ntq),ekc(nstate),freq_r(nw_i:nw),esmr,wtkp,ef
      complex(8) :: zsec(ntq,ntq),zmel(ngb,ntqx,nstate)
      complex(8) :: zw3(ngb,ngb,nwxi:nwx)
C ... Local parameters
      real(8),parameter:: wfaccut=1d-10
      integer itp,itini,itend,iii,iir,it,iw,iws,iwp,ig1,ig2,itpp
      double precision esmrx,wfac,we,omg,wfacx2,weavx2,delwc
      complex(8) :: zz2,zwz3(3),alagr3z2,alagr3z,zdotu,zwkl
      complex(8),allocatable:: zmel1(:),zmel2(:,:),zwk(:,:),zwz5(:,:)
C     complex(8),allocatable:: zwz4(:,:)
C ... for OPENMP
      integer nthreads,omppid,tid,k
!       integer,allocatable:: first0(:),last0(:),idx0(:,:)
!       integer,allocatable:: first5(:),last5(:),idx5(:,:)

      complex(8), parameter :: cmplx0 = 0, cmplx1 = 1
      complex(8), allocatable :: zz(:,:)

      procedure(complex(8)) :: zdotc

      integer(8) :: cr,cm,c1,c2,c10,c20
!      integer :: iwx, wcnt(-nw:nw), wprs(2,256,-nw:nw)

      print *, 'ntq, ntqx:', ntq, ntqx
!       if (ntq /= ntqx) stop 'mksecp: ntq /= ntqx'

      call system_clock(count_rate = cr, count_max = cm)
!      wcnt = 0

! C ... OMP setup
! C$OMP parallel shared(nthreads)
!       nthreads = omppid(0)
! C$OMP END PARALLEL
!       nthreads = 1
!       nthreads = max(nthreads,1)
!       allocate(first0(0:nthreads-1),last0(0:nthreads-1),idx0(2,3*ngb))
!       allocate(first5(0:nthreads-1),last5(0:nthreads-1),idx5(2,3*ntqx))
!       call thdist2(0,0d0,0d0,nthreads,1,3,1,ngb,first0,last0,idx0)
!       call thdist2(0,0d0,0d0,nthreads,1,3,1,ntqx,first5,last5,idx5)



      if (iSigma == 3) then
        allocate(zwz5(3,ntqx))
!         allocate(zwz5(ntqx,3))
        allocate(zwk(ngb,3))
C       allocate(zmel2(ntqx,ngb),zwz5(ntqx,3),zwk(ngb,3))
C       allocate(zmel2(ntqx,ngb),zwz4(3,ntqx),zwk(ngb,3))
        allocate(zmel1(ngb))
        allocate(zz(ntqx,ntqx))
      endif



!       write(*,'("mksecp c2: ",i0)') ntqx

      call system_clock(c10)

!$omp parallel do collapse(2)
      do  itpp = 1, ntqx
        do  itp = 1, ntqx
            zz(itp,itpp) = cmplx0
        end do
      end do
!$omp end parallel do




C --- Big loop over states (q-k,n) ---
C     print *, 'Sigma_en starting time',delwc()
      do  itp = 1, ntqx

C   ... it loop sums contributions from states 1 ... nband
        omg = omega(itp)
        if (omg >= ef) then
          itini = nt0m+1
          itend = nt_max
          iii = 1
        else
          itini = 1
          itend = nt0p
          iii = -1
        endif

        do  it = itini, itend   ! nt0p corresponds to efp

C     ... Make wfac,we
          esmrx = esmr
          if (it<=nctot) esmrx = 0d0
          wfac = wfacx2(omg,ef,ekc(it),esmrx)
          if (GaussSmear()) then
            if (wfac<wfaccut) cycle
            if (it <= nctot .and. wfac/=0) call rx('mksecp: but wfac/=0 for a core state') !faleev
            we = 0.5d0*abs(omg-weavx2(omg,ef,ekc(it),esmr))
          else
            if (wfac == 0d0) cycle
            if (omg >= ef) we = 0.5d0*abs(max(omg-ekc(it),0d0))
            if (omg <  ef) we = 0.5d0*abs(min(omg-ekc(it),0d0))
            if (it <= nctot .and. wfac/=0) call rx('mksecp: but wfac/=0 for a core state') !faleev
          endif
          wfac = iii*wfac*wtkp

C     ... Make iws
          do  iw = 1, nw        !sf 23May02
            iws = iw
            if (freq_r(iw)>we) exit
          enddo


          if (nw_i == 0) then
            if (iws+1>nwx) call rx(' mksecp: energy bin outside range')
          else
            if (omg>=ef .and. iws+1>nwx) then
              write(*,*)'iws+1 nwx=',iws+1,nwx
              call rx(' mksecp: energy bin outside range')
            endif
            if (omg<ef .and. abs(iws+1)>abs(nwxi)) then
              write(*,*)'iws+1 nwxi=',iws+1,nwxi
              call rx(' mksecp: energy bin outside range')
            endif
          endif



C     ... Index iir = +/- 1 depending on whether occ or unocc.  For zw3
          iir = 1
          if (omg < ef .and. nw_i/=0) iir = -1

! Attempt packing as in mkzxq and meltpb, probably has a potential but it is broken for now.
!          do iwp = 1, 3
!            iwx = (iws+iwp-2)*iir
!            if (iwx < -nw .or. iwx > nw) call rx('iwx out of range')
!            wcnt(iwx) = wcnt(iwx)+1
!            if (wcnt(iwx) > 256) call rx('wcnt(iwx) > 256')
!            wprs(1:2,wcnt(iwx),iwx) = [itp,it]
!          end do



!$omp parallel do private(itpp)
            do ig1 = 1, ngb
              zmel1(ig1) = conjg(zmel(ig1,itp,it))
            end do
!$omp end parallel do
            if (iSigma == 3) then

!             do iwp = 1, 3
!                 iw = iws+iwp-2
! !                 call zgemv('t',ngb,ngb,cmplx1,zw3(1,1,iir*iw),ngb,zmel1,1,cmplx0,zwk(1,iwp),1)
!                 call zgemv('t',ngb,ngb,cmplx1,zw3(1,1,iir*iw),ngb,zmel1,1,cmplx0,zwk,1)
!                 call zgemv('t',ngb,ntqx,cmplx1,zmel(1,1,it),ngb,zwk,1,cmplx0,zwz5(iwp,1),3)
!             end do

            iw = min(iir*(iws-1),iir*(iws+1))
!             print *, 'itp,it,iw',itp,it,iw
!             call zgemm('c','n',1,ngb*3,ngb,cmplx1,zmel(1,itp,it),ngb,zw3(1,1,iw),ngb,cmplx0,zwk,1)
            call zgemv('t',ngb,3*ngb,cmplx1,zw3(1,1,iw),ngb,zmel1,1,cmplx0,zwk,1)
            call zgemm('t','n',3,ntqx,ngb,cmplx1,zwk,ngb,zmel(1,1,it),ngb,cmplx0,zwz5,3)

! Swap zwz(iwp) to bring in line with the commented out purely zgemv version. Use only with the zgemm version above.
            if (iir == -1) then
!$omp parallel do private(itpp,zz2)
              do itpp = 1, ntqx
                zz2 = zwz5(1,itpp)
                zwz5(1,itpp) = zwz5(3,itpp)
                zwz5(3,itpp) = zz2
              end do
!$omp end parallel do
            end if

C       ... Calculate contribution to zsec

            if (npm==1) then
!$omp parallel do private(itpp)
              do  itpp = 1, ntqx
                zz(itpp,itp) = zz(itpp,itp) + wfac*alagr3z2(we,freq_r(iws-1),zwz5(:,itpp),itp,itpp)
!                 zsec(itp,itpp) = zsec(itp,itpp) + wfac*alagr3z2(we,freq_r(iws-1),zwz5(:,itpp),itp,itpp)
!                 zsec(itp,itpp) = zsec(itp,itpp) + wfac*alagr3z2(we,freq_r(iws-1),zwz5(itpp,:),itp,itpp)
              enddo
!$omp end parallel do
            else
!$omp parallel do private(itpp)
              do  itpp = 1, ntqx
                zz(itpp,itp) = zz(itpp,itp) + wfac*alagr3z(we,freq_r(iws-1),zwz5(:,itpp))
!                 zsec(itp,itpp) = zsec(itp,itpp) + wfac*alagr3z(we,freq_r(iws-1),zwz5(:,itpp))
!                 zsec(itp,itpp) = zsec(itp,itpp) + wfac*alagr3z(we,freq_r(iws-1),zwz5(itpp,:))
              enddo
!$omp end parallel do
            endif

C     ... QSGW potential iSigma =1,5 (not vectorized)
          else if (iSigma == 1 .or. iSigma==5) then
            zwz3 = cmplx0
            do  iwp = 1, 3
              iw = iws+iwp-2
C             Requires ngb**2 * nt**2 * 3 operations
              do  ig2 = 2, ngb
                zz2 = sum(zmel1(1:ig2-1)*zw3(1:ig2-1,ig2,iir*iw)) + 0.5d0*zmel1(ig2)*zw3(ig2,ig2,iir*iw)
                zwz3(iwp) = zwz3(iwp) + zz2*zmel(ig2,itp,it)
              enddo             !ig2
              zwz3(iwp) = 2d0*dble(zwz3(iwp)) + zmel1(1)*zw3(1,1,iir*iw)*zmel(1,itp,it)
            enddo               !iw
C     ... QSGW potential iSigma =3
C         This contribution to zsec_nn is hermitian
          endif                 ! iSigma - specific evaluation of zsec
        enddo                   ! it

        if (mod(itp,100)==0 .or. itp == ntqx) then
            call system_clock(c20)
            write(*,'("mksecp c20: ",i0,x,g20.12,x,g20.12)')
     &        itp, real(c20 - c10,8)/real(cr,8), real(c20 - c10,8)/real(cr*itp,8)
         end if
      enddo                     !itp
C     print *, 'Sigma_en ending time',ngb,ntqx,delwc()

!$omp parallel
!$omp do private(itpp,itp)
      do  itpp = 1, ntqx
        do  itp = 1, ntqx
            zsec(itp,itpp) = zsec(itp,itpp) + zz(itpp,itp)
        end do
      end do
!$omp end do
!$omp end parallel


      if (iSigma == 3) then
        deallocate(zwz5,zwk)
        deallocate(zmel1)
        deallocate(zz)
C       deallocate(zmel2,zwz5,zwk)
C       deallocate(zmel2,zwz4,zwk)
      endif

!       deallocate(first0,last0,idx0)
!       deallocate(first5,last5,idx5)

C      call zprm('zsec',zsec,ntq,ntqx,ntqx)
C      stop


!       do iw = -nw,nw
!          if (wcnt(iw) == 0) cycle
!          print *, 'iw, wcnt(iw), itps,it', iw, wcnt(iw)
!          do it = 1, wcnt(iw)
!             print *, wprs(1:2,it,iw)
!          end do
!       end do

      end

C     Attempt at level 3 blas vectorization
C      subroutine mksecp(iSigma,npm,nstate,ntq,ntqx,nt0p,nt0m,nt_max,
C     .  nctot,ngb,ekc,ef,omega,GaussSmear,esmr,nw,wtkp,nwxi,nwx,nw_i,
C     .  freq_r,zw3,zmel,zsec)
CC- Adds pole contribution to self energy for QSGW Vxc
CC ----------------------------------------------------------------------
CCi Inputs
CCi   iSigma: mode determining form of static approximation to SE
CCi         : 0 SE_nn'(ef)+imag integral:delta_nn'(SE_nn(e_n)-SE_nn(ef))
CCi         : 1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
CCi         : 2 SE_nn'((e_n+e_n')/2)
CCi         : 3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
CCi         : 4 SE_nn'(ef)
CCi         : 5 delta_nn' SE_nn(e_n)
CCi   npm   :1 Use time reversal symmetry
CCi         :2 No time reversal symmetry
CCi   nstate:number of states contributing to sigma
CCi   ntq   :dimension of 1-particle hamiltonian; number of eigenstates t
CCi         :dimensions zsec
CCi   ntqx  :Generate SE for ntqx states
CCi   nt0p  :upper bound to it index when omg(itp)<ef
CCi   nt0m  :lower bound to it index when omg(itp)>ef
CCi   nt_max:upper bound to it index when omg(itp)>ef
CCi   nctot :number of core states
CCi   ngb   :Number of product basis functions
CCi   ekc   :ekc(1:nstate) QP levels
CCi   ef    :Fermi level
CCi   omega :frequencies for each QP level
CCi   GaussSmear
CCi         : T, Broaden poles in G with Gaussian of width esmr
CCi         : F, Broaden poles in G with rectangle of width esmr
CCi   esmr  :Broadening of G in integration; see wintzg
CCi   nw    :Number of imaginary mesh points
CCi   nw_i  :lower bound of freq_r
CCi   nwxi  :dimensions zw3 - zw3(:,:,nwxi:nwx)
CCi   nwx   :dimensions zw3 - zw3(:,:,nwxi:nwx)
CCi   wtkp  :k-point weight
CCi   freq_r:Frequencies on the real axis
CCi   zmel  :matrix elments of the product basis
CCi   zw3   :
CCo Outputs
CCo   zsec  :Pole contribution to S_ij= <i|Re[S](e_i)|j> for a given kp
CCs Command-line switches
CCl Local variables
CCl   iir   : 1 if current omega is > ef
CCl         :-1 if current omega is < ef
CCr Remarks
CCr  This is a vectorized version of mksecp0.
CCr  Implemented for iSigma=3 only so far.
CCr
CCr  Extra information is required for frequencies not needed.
CCr  Not worth the extra effort
CCu Updates
CCu   22 May 13 Adapted from Takao's sxcf_fal2.sc.F
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      logical GaussSmear
C      integer iSigma,npm,nstate,ntq,ntqx,nt0m,nt0p,nt_max,nctot,ngb,
C     .        nw,nwxi,nwx,nw_i
C      real(8) :: omega(ntq),ekc(nstate),freq_r(nw_i:nw),esmr,wtkp,ef
C      complex(8) :: zsec(ntq,ntq),zmel(ngb,nstate,ntqx)
C      complex(8) :: zw3(ngb,ngb,nwxi:nwx)
CC ... Local parameters
C      real(8),parameter:: wfaccut=1d-10
C      integer itp,itini,itend,iii,it,iw,ixs,iir,ix0,ix,ig1,ig2,itpp
C      double precision esmrx,omg,
C     .  wfacx2,weavx2 !,delwc
C      complex(8) :: zz2,zwz3(3),alagr3z2,alagr3z
C      complex(8),allocatable:: zmel1(:),zwz4(:,:),zwk(:),zwk2(:,:,:)
C      real(8),allocatable:: wfac(:,:),we(:,:)
C      integer,allocatable:: iwi(:,:),iti(:,:)
C      integer iwl(3*nstate,2),iwp,nwp,nwl,ntx
C
C      if (iSigma /= 3) call rx('mksecp implemented only for isigma=3')
C
C      allocate(zmel1(ngb))
CC     if (iSigma == 3) allocate(zwz4(3,ntqx),zwk(ngb))
C      if (iSigma == 3) allocate(zwz4(3,ntqx),zwk2(ngb,ntqx,3))
C
C --- First pass through loop over states (q-k,t') ---
C     This pass finds all (t,t') pairs that contribute to zsec
C     Also makes wfac(it,itp),we(it,itp),ixs(it,itp)
C      allocate(iwi(nw,ntqx),iti(nstate,ntqx))
C      allocate(wfac(nstate,ntqx),we(nstate,ntqx))
C      iwi = -1 !loop sets iwi(iw)=1 if central bin needed for it,itp
C      call iinit(iti,nstate*ntqx) ! iti(t,t')=1 if pair contributes to zsec
C      do  itp = 1, ntqx
CC       print *, 'Sigma_en starting time',itp,delwc()
C        omg = omega(itp)
C        if (omg >= ef) then
C          itini = nt0m+1
C          itend = nt_max
C          iii = 1
C        else
C          itini = 1
C          itend = nt0p
C          iii = -1
C        endif
CC       Index iir for zw3
C        iir = 1
C        if (omg < ef .and. nw_i/=0) iir = -1
C
CC   ... Inner loop t'
C        do  it = itini, itend
C
CC     ... Make wfac,we
C          esmrx = esmr
C          if (it <= nctot) esmrx = 0d0
C          wfac(it,itp) = wfacx2(omg,ef,ekc(it),esmrx)
C          if (GaussSmear()) then
C            if (wfac(it,itp)<wfaccut) cycle
C            if (it <= nctot .and. wfac(it,itp)/=0) then !faleev
C              call rx('mksecp: but wfac/=0 for a core state')
C            endif
C            we(it,itp) = 0.5d0*abs(omg - weavx2(omg,ef,ekc(it),esmr))
C          else
C            if (wfac(it,itp) == 0d0) cycle
C            if (omg >= ef) we(it,itp) = 0.5d0*abs(max(omg-ekc(it),0d0))
C            if (omg <  ef) we(it,itp) = 0.5d0*abs(min(omg-ekc(it),0d0))
C            if (it <= nctot .and. wfac(it,itp)/=0) then !faleev
C              call rx('mksecp: but wfac/=0 for a core state')
C            endif
C          endif
C          iti(it,itp) = 1
C          wfac(it,itp) = iii*wfac(it,itp)*wtkp
C
CC     ... Make iwi
C          do  iw = 1, nw        !sf 23May02
C            iwi(it,itp) = iw
C            if (freq_r(iw)>we(it,itp)) exit
C          enddo
C          if (nw_i == 0) then
C            if (iw+1>nwx)
C     .      call rx(' mksecp: energy bin outside range')
C          else
C            if (omg >= ef .and. iwi(it,itp)+1> nwx) then
C              write(*,*)'iw+1 nwx=',iwi(it,itp)+1,nwx
C              call rx(' mksecp: energy bin outside range')
C            endif
C            if (omg < ef .and. abs(iwi(it,itp)+1)> abs(nwxi) ) then
C              write(*,*)'iw+1 nwxi=',iwi(it,itp)+1,nwxi
C              call rx(' mksecp: energy bin outside range')
C            endif
C          endif
C        enddo
C      enddo
C
CC --- Second pass over states (q-k,t') evaluates zsec  ---
C      do  itp = 1, ntqx
C
C        omg = omega(itp)
C        if (omg >= ef) then
C          itini = nt0m+1
C          itend = nt_max
C          iii = 1
C        else
C          itini = 1
C          itend = nt0p
C          iii = -1
C        endif
C
CC   ... Make iwl = list of all energies in zw3 required for this (tt') pair
CC       nwl = length(iwl)
C        iwp = 0
C        do  it = itini, itend
C          if (iti(it,itp) == 0) cycle
C          do  ix0 = 1, 3
C            ix = iwi(it,itp)+ix0-2
C            iwp = iwp+1
C            iwl(iwp,2) = iir*ix
C          enddo
C        enddo
C        nwp = iwp
C        call ivheap(1,nwp,iwl(1,2),iwl(1,1),0)
C        iwl(1,1) = iwl(1,2)
C        nwl = 1
C        do  iwp = 1, nwp
C          if (iwl(iwp,2) /= iwl(nwl,1)) then
C            nwl = nwl+1
C            iwl(nwl,1) = iwl(iwp,2)
C          endif
C        enddo
C
CC   ... Contracted part of zmel
C        ntx = sum(iti(itini:itend,itp))
C        stop
CC       allocate(zm(ngb,ntx))
CC        ix = 0
CC        do  it = itini, itend
CC          if (iti(it,itp) == 0) cycle
CC          ix = ix+1
CC          zm(:,ix) = zmel(:,it,itp)
CC        enddo
CC        if (ntx /= ix) call rx('bug in mksecp')
C
CC   ... Make cc(t')_g,t,iw
CC        allocate(cc(ntx,ngb,nwl))
CC        do  iw = 1, nwl
CCC          iwp = iwl(iw,1)
CCC          call zcopy(ngb*ngb,zw3(1,1,iwp),1,zw4(1,1,iw))
CC          call zgemm('C','N',ntx,ngb,ngb,(1d0,0d0),zm,ngb,
CC     .      zw3(1,1,iwp),ngb,(0d0,0d0),cc(1,1,iw),ntx)
CC        enddo
C
C
CC             Calculate zwz4(ix0,tpp,t,tp)
CC             zwz4 = sum_t,tp,g1,g2,ix zm+(g1,t,tp) zw3(g1,g2,iir*ix) zm(g2,t,tpp)
CC                  = sum_t,tp,g2 zwk(g2)*zmel(g2,t,tpp)
CC             Requires ngb * nt**3 * 3 operations
C        do  ix0 = 1, 3
C          ix = ixs+ix0-2
C          do  itpp = 1, ntqx
C            do  ig2 = 1, ngb
C              zwz4(ix0,itpp) = zwz4(ix0,itpp) +
C     .          zwk(ig2)*zmel(ig2,it,itpp)
C                enddo           !ig2
C              enddo             !itpp
C            enddo               !ix0
C
C
C
C        print *, itp,iw
CC       print *, iwl(1:iw,1)
C
C      enddo
C
C      stop
C
C
C
C      do  itp = 1, ntqx
CC       print *, itp,sum(iti(:,itp))
C        print '(i3,2x,50i3)', itp,iwi(:,itp)
C      enddo
C      stop
C
C      end
