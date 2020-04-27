      subroutine hstra(MOL,pv,ldip,strx,drstrx,nlmf,nlm,nlmq1,nlmq,hl,
     .                 cg,indxcg,jcg,vol)
C- Make structure constants from reduced strux at energy zero taking
C- advantage of symmetry properties of zero energy strux
C ----------------------------------------------------------------------
Ci Inputs:
Ci  MOL   : if T skip dipole correction (use non-periodic setup)
Ci  pv    : if T calculate pressure
Ci  ldip  : 3 include 'spherical' dipole correction to Ewald
Ci        : 2 include 'slab' dipole correction to Ewald
Ci        : any other number - skip the dipole correction
Ci  nlmf,nlm :make B_LL' for L = 1, nlmf, L'= 1, nlm
Ci  nlmq1 : leading dimensions of B, nlmq1=(ll(nlmq)+1)**2
Ci  nlmq  : max L-cutoff for multipoles, leading dimension of B'
Ci  hl    : solid Hankels for given point tau
Ci  cg,indxcg,jcg : Clebsch Gordan coefficients in condensed form
Ci  vol   : unit cell volume
Co Outputs:
Co  strx  : coefficients B_LL'(tau) (structure constants)
Co  dstrx : derivatives of structure constants (x dB_LL'/dx) at x = tau
Co          if pv = F dstrx is not touched
Cr Remarks:
Cr  strx are zero energy structure constants B_LL' are calculated as:
Cr    B_LL'(tau) = 4\pi \sum_L" (-1)^l' C_LL'L" H_l"(tau)
Cr  where C_LL'L" are Gaunt coefficients (C_LL'L" = \int Y_L * Y_L' * Y_L")
Cr  and H_L are solid Hankels at zero energy
Cr    H_L(tau) = h_l"(|tau|) Y_L"(tau)
Cr  prepared by either rcnsl0.f (periodic branch, MOL = .F.) or
Cr  soldh.f (non-periodic branch, MOL = .T.).
Cr
Cr  If pv is set, returns tau*d/dtau B in drstrx
Cr
Cr  The program is an accelerated version of hstr.f which takes
Cr  advantage of the symmetry of B_LL':
Cr          B_L'L(tau) = B_LL'(tau) * (-1)^(l+l')   (*)
Cr  (*), in turn, is a consequence of the property
Cr  Y_L(-tau) = (-1)^l Y_L(tau)
Cr
Cr  hstra first calculates the lower left triangle of B_LL'(L = 1, nlmf,
Cr  L' = 1, L) and then fills up the upper right triangle using (*).
Cr  For the sake of clarity, the program is restricted to nlmf >= nlm.
Cr  No such restriction is imposed in hstr.f though.
Cr
Cu Updates
Cu    05 Mar 2010 (SL) created from hstr.f
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      logical, intent(in) :: MOL,pv
      integer, intent(in) :: nlmq1,nlmq,nlm,nlmf,indxcg(*),jcg(*),ldip
      double precision, intent(in) :: cg(*),hl(*),vol
      double precision, intent(out) :: strx(nlmq1,nlmq),
     .                                 drstrx(nlmq,nlmq)
C Local Parameters
      integer lmxx,icg,icg1,icg2,ii,ilm,indx,klm,l,lk,ll,llm,lm,lp,
     .  lmax,mlm,klmx
      parameter (lmxx=12)
      double precision sig(0:lmxx),fourpi,fpibv,sum,sumr

!       call tcn('hstra: make strux')

      fourpi = 16d0*datan(1d0)
      lmax = ll(nlmf) + ll(nlm)
      if (lmax > lmxx) call rx0(' change dimensions in hstra')
      if (nlm > nlmf) call rx0(' hstra: only works with '//
     .  'nlmf >= nlm. Use hstr instead.')

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
        klmx = min0(mlm,nlm)
        do  klm = 1, klmx
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
C...        Only lp = lm + lk contribute
            if (lp /= lm+lk) cycle
            sum = sum + cg(icg)*hl(llm)
            if (pv) sumr = sumr - (lp+1) * cg(icg)*hl(llm)
          enddo
          strx(mlm,klm) = sum*fourpi*sig(lk)
          if (pv) then
C...        need to cut off l+1 terms for B'. Quick fix for now
            if (mlm <= nlmq) drstrx(mlm,klm) = sumr*fourpi*sig(lk)
          endif
        enddo
      enddo
C --- Add the remaining off-diagonal terms using B_LL' = -1^(l+l')*B_L'L
      if (nlm > 1) then
        do  mlm = 1, nlm-1
          lm = ll(mlm)
          do  klm = mlm+1, nlm
            lk = ll(klm)
            strx(mlm,klm) =  (sig(lm)*sig(lk))*strx(klm,mlm)
            if (pv) drstrx(mlm,klm) = (sig(lm)*sig(lk))*drstrx(klm,mlm)
          enddo
        enddo
      endif

C --- the following includes extra p terms 'implicitly' ---
      if ((ldip == 2 .or. ldip == 3) .and. (.not. MOL)) then
        if (nlm > 1) then
          fpibv = fourpi/vol
          do  ilm = 2, 4
            strx(ilm,ilm) = strx(ilm,ilm) - fpibv
          enddo
          if (pv) then
            do  ilm = 2, 4
               drstrx(ilm,ilm) = drstrx(ilm,ilm) + 3d0*fpibv
            enddo
          end if
        endif
      endif

      call tbshfl(0,nlmq1,nlmf,nlm,strx)
      call strfac(0,nlmq1,nlmf,nlm,strx)

      if (pv) then
         call tbshfl(0,nlmq,nlm,nlm,drstrx)
         call strfac(0,nlmq,nlm,nlm,drstrx)
      end if

!       call tcx('hstra: make strux')
      end subroutine hstra
