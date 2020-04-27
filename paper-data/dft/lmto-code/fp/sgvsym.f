      subroutine sgvsym(ngrp,g,ag,ng,gv,ips0,bgv)
C- Setup for symmetrization of a function in Fourier representation.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ngrp    number of group operations
Ci           ngrp=-2 => special AFM symmetrization
Ci           Only 2nd group operation contributes to bgv.
Ci   g,ag    space group
Ci   ng,gv   list of reciprocal lattice vectors
Co Outputs:
Co   ips0    pointer to first vector in star of this vector
Co   bgv     phase factor sum; see Remarks
Cr Remarks:
Cr   The reciprocal lattice vectors are assumed to be sorted by length
Cr   Adapted from nfp su_gvsym.f
Cr
Cr   Symmetrized f means  f(G) = f(g(G)).
Cr   Let G be all r.l.v and G' be just first members of each star.  Then
Cr
Cr   f(r) = sum_G f(G) exp(i G.r)
Cr        = sum_G' sum_(G in star of G')  f(G') exp(i G.r)
Cr        = sum_G' f(G') exp(i G.r) sum_(G in star of G') exp(i(G-G').r)
Cr        = sum_G' f(G') exp(i G.r) bgv(star)
Cb Bugs
Cb   This routine copies su_gvsym.f conventions, which use inv(g)
Cu Updates
Cu    7 Sep 98 Adapted from nfp su_gvsym.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ngrp,ng,ips0(ng)
      double precision g(3,3,ng),ag(3,ng),gv(ng,3)
      double complex bgv(ng)
C ... Local parameters
      integer i,i00,irep,i0,nstar,k,j,j0,iprint,ksum,kstar,stdo,lgunit
      double precision tpi,df,scalp,gg0,gg,fac,vv,v(3)

      stdo = lgunit(1)
      tpi = 8d0*datan(1d0)

      do  i = 1, ng
        ips0(i) = 0
        bgv(i) = (0d0,0d0)
      enddo

C --- Loop over stars of k ---
      i00 = 1
      do  irep = 1, ng+1
        i0 = 0
        do  i = i00, ng         ! Look for next unclassified vector
          i0 = i
          if (ips0(i) == 0) goto 5
        enddo
        goto 10                 ! No more vectors
    5   continue

C   --- Apply all point ops, find in list, add to phase sum ---
        nstar = irep
        do  k = 1, abs(ngrp)
C       do  k = 1, max(1,ngrp)
C     ... Find G' = g(k) G; j0 is index to G'
          v(1) = g(1,1,k)*gv(i0,1)+g(1,2,k)*gv(i0,2)+g(1,3,k)*gv(i0,3)
          v(2) = g(2,1,k)*gv(i0,1)+g(2,2,k)*gv(i0,2)+g(2,3,k)*gv(i0,3)
          v(3) = g(3,1,k)*gv(i0,1)+g(3,2,k)*gv(i0,2)+g(3,3,k)*gv(i0,3)
          do  j = i0, ng
            df = (v(1)-gv(j,1))**2+(v(2)-gv(j,2))**2+(v(3)-gv(j,3))**2
            j0 = j
            if (df < 1d-8) goto 6
          enddo
          write (stdo,1) i0,k,gv(i0,1),gv(i0,2),gv(i0,3),v
    1     format(' ---- vec',i6,'   op',i3,2x,3F8.4,'  ->',3F8.4)
          call rxi('SGVSYM: cannot find mapped vector in list:',i0)
    6     continue
          ips0(j0) = i0
          if (k == 1 .and. ngrp < 0) cycle
          scalp = gv(j0,1)*ag(1,k)+gv(j0,2)*ag(2,k)+gv(j0,3)*ag(3,k)

          if (ngrp > 0) then
            bgv(j0) = bgv(j0) + cdexp(dcmplx(0d0,tpi*scalp))
          else
            bgv(j0) = cdexp(dcmplx(0d0,tpi*scalp))
            scalp = gv(i0,1)*ag(1,k)+gv(i0,2)*ag(2,k)+gv(i0,3)*ag(3,k)
            bgv(i0) = cdexp(dcmplx(0d0,tpi*scalp))
          endif
        enddo
        i00 = i0
      enddo
      call rxi('SGVSYM: this cannot happen. irep=',irep)
   10 continue

      if (iprint() >= 20) call awrit2(' SGVSYM: %i symmetry stars'//
     .  ' found for %i reciprocal lattice vectors',' ',80,stdo,nstar,ng)

C --- Multiply phase sums by (star order)/(group order) ---
      ksum = 0
      do  i0 = 1, ng
        if (ips0(i0) == i0) then ! Block for a star: condition => 1st vector
          kstar = 0 ! kstar will be size of star
          gg0 = gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2
C         Count number of elements in star
          do  i = i0, ng
            if (ips0(i) == i0) kstar = kstar+1
            gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
            if (dabs(gg-gg0) > 1d-3) exit
          enddo
          ksum = ksum+kstar ! Total number of elements
C         Scale by the number in star
          if (ngrp == -2) then ! Special AFM symmetrization
            fac = 1
          else
            fac = dble(kstar)/dble(ngrp)
          endif
          do  i = i0, ng
            if (ips0(i) == i0) bgv(i) = bgv(i)*fac
            gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
            if (dabs(gg-gg0) > 1d-3) exit
          enddo
        endif
      enddo
      if (ksum /= ng) call rxi('SGVSYM error, ksum=',ksum)

C      if (ngrp < 0) then
C        do  i = 1, ng
C          j = ips0(i)
C          if (j == i) cycle
CC          print *, i, sngl(gv(i,:))
CC          print *, j, sngl(gv(j,:))
C          print *, i,j, sum(abs(sngl(gv(j,:)+gv(i,:))))
C        enddo
C        stop
C      endif

C --- Printout ---
      if (iprint() < 60) return
      j = min(ng,300)
      if (iprint() >= 90) j = ng
      nstar = 0
      write (stdo,2)
    2 format(//'   elt     no           vector                phase')
      do  i0 = 1, j
        if (ips0(i0) == i0) then
          nstar = nstar+1
          vv = dsqrt(gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2)
          write (stdo,3) nstar,vv
    3     format(/'   Star',i6,'  length =',f8.4)
          kstar = 0

          do  i = i0, ng
            if (ips0(i) == i0) then
              kstar = kstar+1
              write (stdo,4) kstar,i,gv(i,1),gv(i,2),gv(i,3),bgv(i)
    4         format(i5,i8,2x,3F7.2,3x,2F7.2)
            endif
          enddo
        endif
      enddo

      end
