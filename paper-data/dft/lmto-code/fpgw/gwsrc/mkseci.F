      subroutine mkseci(iw0,iSigma,npm,nstate,ntq,ntqx,nctot,ef,omega,
     .  ekc,GaussSmear,esmr,nw,freqx,wx,a,expa,wtkp,zwz,zwz0,zsec)
C- Makes or adds part of frequency integration for sigma
C ----------------------------------------------------------------------
Ci Inputs
Ci   iw0   : 0 Perform omega-integration for entire quadrature of nw points.
Ci         :   In this case zwz is dimensioned zwz(:,:,:,1:nw)
Ci         : 1 Add the contribution to omega-integration for point ix0
Ci         :   In this case zwz is dimensioned zwz(:,:,:,1)
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
Ci   ef    :Fermi level
Ci   omega :frequencies for each QP level
Ci   ekc   :ekc(1:nstate) QP levels
Ci   GaussSmear
Ci         : T, Broaden poles in G with Gaussian of width esmr
Ci         : F, Broaden poles in G with rectangle of width esmr
Ci   esmr  :Broadening of G in integration; see wintzg
Ci   nw    :Number of imaginary mesh points
Ci   freqx :Vector of frequencies transformed to (0,1); see wintz0 or wintzg
Ci   wx    :Integration weights at gaussian points freqx between (0,1)
Ci   a     :Inverse width of gaussian to be integrated out analytically,
Ci         :making residual small; see wintzg
Ci   expa  :Vector of exp(-a^2 w'^2) for the nw frequencies iw'
Ci   wtkp  :k-point weight
Ci   zwz0  :Matrix elements of static correlation potential, zmel*(W(0)-v)*zmel
Ci   zwz   :Matrix elements of dynamic correlation potential, zmel*(W(iw)-v)*zmel
Ci         :If iw0=0, zwz = zwz(:,:,:,nw) --- contains zwz for all energies
Ci         :If iw0>0, zwz = zwz(:,:,:,1) --- contains zwz for iw0 only
Co Outputs
Co   zsec  :S_ij= <i|Re[S](e_i)|j>
Cs Command-line switches
Cl Local variables
Cl    iw   :index to current energy point
Cl    ix   :index to position in zwz containing zmel*(W(iw)-v)*zmel
Cr Remarks
Cr   This routine calculates <i|Re[S](e_i)|j>
Cr   Im-axis integral should give Hermitian part of this integral.
Cr   But because of e_i there is a difference between
Cr   <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
Cr   Prior to Jan2006 sxcf_fal2.sc.F symmetrized it explicity
Cr   but is not done here now because it is done hqpe.sc.F anyway
Cu Updates
Cu   10 May 13  Adapted from Takao's sxcf_fal2.sc.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical GaussSmear
      integer iw0,iSigma,nstate,ntq,ntqx,nw,nctot,npm
      double precision esmr,a,wtkp,ef
      real(8) :: omega(ntqx),ekc(nstate),freqx(nw),wx(nw),expa(nw)
      complex(8) :: zwz(ntqx,ntqx,nstate,nw),zwz0(ntqx,ntqx,nstate)
      complex(8) :: zsec(ntq,ntq)
C ... Local parameters
      integer itp,itpp,it,iw,ix,low,high
      double precision omega0,we,esmrx
      complex(8):: img=(0d0,1d0),zwzs0,zs,wintzg,wintzav
      complex(8),allocatable:: zwzs(:,:),zseci(:)

!       integer(8) :: cr,cm,c1,c2

!       call system_clock(count_rate = cr, count_max = cm)
!       call system_clock(c1)
!       print *,'mkseci: isigma, gausssmear, iw0', isigma, gausssmear, iw0
C ... Setup
      if (iw0 == 0) then
        low = 1; high = nw
      else
        low = iw0; high = iw0
      endif

C --- Diagonal part, applies to iSigma = 0,1,5 ---
      if (iSigma==0 .or. iSigma==1 .or. iSigma==5) then
        allocate(zwzs(nw,2),zseci(nstate))
        do  itp = 1, ntqx
          do  it = 1, nstate
            esmrx = esmr; if (it <= nctot) esmrx = 0d0
            we = (omega(itp) - ekc(it))/2
            if (GaussSmear) then
              do  iw = low, high
                if (iw0 == 0) then
                  ix = iw
                else
                  ix = 1
                endif
                zwzs(iw,1) = dble(zwz(itp,itp,it,ix)) ! w(iw) + w(-iw) symmetric part
                if (npm==2) then
                  zwzs(iw,2) = dimag(zwz(itp,itp,it,ix)) ! w(iw) - w(-iw)
                endif
              enddo
C             zwzi(it,itp,itp) =
              zseci(it) = wintzg(npm,iw0,zwzs,zwz0(itp,itp,it),freqx,wx,a,expa,we,nw,esmrx)
            else
              call rx('mkseci: GaussSmear off not yet supported')
              if (npm==2 .or. iw0 /= 0)
     .          call rx('mkseci: GaussSmear off implemented only with TR symmetry')
C             zwzi(it,itp,itp) =
              zseci(it) = wintzav(zwz(itp,itp,it,1:nw),zwz0(itp,itp,it),freqx,wx,a,expa,we,nw,esmrx)
            endif
          enddo                 !it

C         Add to Sigma_{ii}(e_i):
C         Sum over occ + unoccupied states,  and multiply by kp weight
          zs = sum(zseci(:))
C         zs = sum(zwzi(:,itp,itp))
C         zs = (zs + dconjg(zs))/2    !Make hermitian
          zsec(itp,itp) = zsec(itp,itp) + wtkp*zs
        enddo                   !itp
        deallocate(zwzs,zseci)
        if (iSigma==5) return     !Nothing more to do in diagonal-only case
      endif                     ! Diagonal part

C --- All iSigma modes other than 5 ---
!$OMP parallel private (itp,itpp,it,iw,ix,omega0,we,esmrx,zwzs0,zs,zwzs,zseci)
      allocate(zwzs(nw,2),zseci(nstate))
!$OMP do
      do  itpp = 1, ntqx
        do  itp = 1, ntqx
          do  it = 1, nstate
            if (iSigma==0 .or. iSigma==1 .or. iSigma==4) then
              omega0 = ef
            elseif (iSigma==2) then
              omega0 = .5d0*(omega(itp)+omega(itpp))
            elseif (iSigma==3) then
              omega0 = omega(itp)
            else
              call rx("mkseci: iSigma must be between 0 and 5")
            endif
            esmrx = esmr; if (it <= nctot) esmrx = 0d0
            we = (omega0 - ekc(it))/2
            if (GaussSmear) then
              do  iw = low, high
                if (iw0 == 0) then
                  ix = iw
                else
                  ix = 1
                endif
C               Symmetric (hermitian) part [w(iw) + w(-iw)]/2
!                 zwzs(iw,1) = (zwz(itp,itpp,it,ix) + dconjg(zwz(itpp,itp,it,ix)))/2
                zwzs(iw,1) = 0.5d0*(zwz(itp,itpp,it,ix) + dconjg(zwz(itpp,itp,it,ix)))
!                 zwzs(iw,1) = real(zwz(it,itp,itpp,ix))
C               Antisymmetric part ! [w(iw) - w(-iw)]/2i
                if (npm == 2) then
! the -1 is gone to compensate for 'c,n' type zgemm in matzwz3
                  zwzs(iw,2) = (0d0,0.5d0)*(zwz(itp,itpp,it,ix) - dconjg(zwz(itpp,itp,it,ix)))
!                   zwzs(iw,2) = (zwz(itp,itpp,it,ix) - dconjg(zwz(itpp,itp,it,ix)))/2/img
!                   zwzs(iw,2) = -(0d0,0.5d0)*(zwz(itp,itpp,it,ix) - dconjg(zwz(itpp,itp,it,ix)))
!                   zwzs(iw,2) = aimag(zwz(it,itp,itpp,ix))
                endif
              enddo
              zwzs0 = (zwz0(itp,itpp,it) + dconjg(zwz0(itpp,itp,it)))/2 !not necessary
C             zwzi(it,itp,itpp) =
              zseci(it) = wintzg(npm,iw0,zwzs,zwzs0,freqx,wx,a,expa,we,nw,esmrx)
            else
              call rx('mkseci: GaussSmear off not yet supported')
              if (npm==2 .or. iw0 /= 0)
     .          call rx('mkseci: GaussSmear off implemented only with TR symmetry')
              zseci(it) =
     .          wintzav(zwz(itp,itpp,it,1:nw),zwz0(itp,itpp,it),freqx,wx,a,expa,we,nw,esmrx)
            endif
          enddo                 ! loop over intermediate states

C         Add to Sigma_{ij}(e_i):
C         Sum over occ + unoccupied states, and multiply by kp weight
C         Skip cases diagonal part already done at top of code
          if ((iSigma==1 .or. iSigma==0) .and. itp==itpp) cycle
C         zs = sum(zwzi(:,itp,itpp))
          zs = sum(zseci(:))
C         Prior to 2006
C         zs = (sum(zwzi(:,itp,itpp))+dconjg(sum(zwzi(:,itpp,itp))))/2
          zsec(itp,itpp)  = zsec(itp,itpp) + wtkp*zs
        enddo                   !itp
      enddo                     !itpp
!$OMP end do
      deallocate(zwzs,zseci)
!$OMP end parallel

!       call system_clock(c2)
!       write(*,'("mseci:",x,g20.12)') real(c2 - c1,8)/real(cr,8)


      end
