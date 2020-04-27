      subroutine hstr(MOL,pv,ldip,strx,drstrx,nlmf,nlm,nlmq1,nlmq,hl,cg,
     .                indxcg,jcg,vol)
C- Make structure constants from reduced strux at energy zero
C ----------------------------------------------------------------------
Ci Inputs:
Ci  MOL   : if T skip dipole correction (use non-periodic setup)
Ci  pv    :if T calculate pressure
Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
Ci        : 2 include 'slab' dipole correction to Ewald
Ci        : any other number - skip the dipole correction
Ci  nlmf,nlm :make B_LL' for L = 1, nlmf, L'= 1, nlm
Ci  nlmq1 :leading dimensions of B, nlmq1=(ll(nlmq)+1)**2
Ci  nlmq  :max L-cutoff for multipoles, leading dimension of B'
Ci  hl    :radial Hankels for given point |tau|
Ci  cg,indxcg,jcg : Clebsch Gordan coefficients in condensed form
Ci  vol   :unit cell volume
Co Outputs: strx,drstrx
Co  strx  :coefficients B_LL'(tau) (structure constants)
Co  dstrx :derivatives of structure constants (x dB_LL'/dx) at x = tau
Co         if pv = F dstrx is not touched
Cr Remarks: B_LL'(tau) = 4\pi \sum_L" (-1)^l' H_l"(tau) Y_L"(tau)
Cr         HY are reduced structure constants (see rcnsl0, soldh)
Cr         If pv is set, returns tau*d/dtau B in drstrx
Cu Updates
Cu   28 Apr 10 (SL) Slab dipole correction to Ewald
C ----------------------------------------------------------------------
      implicit none
      integer nlmq1,nlmq,nlm,nlmf,indxcg(*),jcg(*),ldip
      double precision cg(*),strx(nlmq1,nlm),drstrx(nlmq,nlm),hl(*),vol
      integer lmxx,icg,icg1,icg2,ii,ilm,indx,ipow,klm,l,lk,ll,llm,lm,lp,
     .  lmax,mlm
      parameter (lmxx=12)
      double precision sig(0:lmxx),fourpi,fpibv,sum,sumr
      logical MOL,pv

!       call tcn('hstr: make strux')
      fourpi = 16d0*datan(1d0)
      lmax = ll(nlmf) + ll(nlm)
      if (lmax > lmxx) call rx0(' change dimensions in hstr')
C --- (-1)^l ---
      sig(0) = 1d0
      if (lmax > 0) then
        do  l = 1, lmax
          sig(l) = -sig(l-1)
        enddo
      endif
C --- add together Gaunt-sums ---
      do  mlm = 1, nlmf
        lm = ll(mlm)
        do  klm = 1, nlm
          lk = ll(klm)
          sum = 0d0
          sumr = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2+min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  icg = icg1, icg2
            llm = jcg(icg)
            lp = ll(llm)
            ipow = (lm + lk - lp)/2
            if (ipow == 0) then
              sum = sum + cg(icg)*hl(llm)
              if (pv) then
                sumr = sumr - (lp+1) * cg(icg)*hl(llm)
              endif
            endif
          enddo
          strx(mlm,klm) = sum*fourpi*sig(lk)
          if (pv) then
c need to cut off l+1 terms for B'. Quick fix for now
            if (mlm <= nlmq) drstrx(mlm,klm) = sumr*fourpi*sig(lk)
          endif
        enddo
      enddo
C --- the following includes extra p terms 'implicitly' ---
      if ((ldip == 2 .or. ldip == 3) .and. (.not. MOL)) then
        if (nlm > 1) then
          fpibv = fourpi/vol
          do  ilm = 2, min0(4,nlm,nlmf)
            strx(ilm,ilm) = strx(ilm,ilm) - fpibv
            if (pv) then
              drstrx(ilm,ilm) = drstrx(ilm,ilm) + 3d0*fpibv
            endif
          enddo
        endif
      endif

      call tbshfl(0,nlmq1,nlmf,nlm,strx)
      call strfac(0,nlmq1,nlmf,nlm,strx)

      if (pv) then
         call tbshfl(0,nlmq,nlmq,nlm,drstrx)
         call strfac(0,nlmq,nlmq,nlm,drstrx)
      end if

!       call tcx('hstr: make strux')
      end
