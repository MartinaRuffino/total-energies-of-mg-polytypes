      subroutine rotycs(mode,a,nbas,nsp,lmaxu,s_spec,s_site,lldau)
C- Rotate n LDA+U formatted matrix (densmat or vorb) to/from real, spherical harmonics
C-------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa idu
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   mode  :1  converts a from real to spherical
Ci         :-1 converts a from spherical to real
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Cio Inputs/Outputs
Ci   a     :a(:,:,isp,iblu) is rotated in place
Cr Remarks
Cr This uses the standard convention (real R to spherical Y)
Cr   ( Y_l,m  ) =           ((-1)^m + (-1)^m i) (R_l,m )
Cr   (        ) = 1/sqrt(2) (                 ) (      )
Cr   ( Y_l,-m ) =           (1      - i       ) (R_l,-m)
Cr Or equivalently
Cr Yl,-m = (Rlmc-iRlms)/sqrt(2)  Ylm = (-1)**m*conjg(Yl,-m)
Cr Order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
Cr                             ... Rl,-m ...      ... Rl,m ...
Cr Order of spherical harmonics m = -l:l
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   18 Jan 06 A. Chantis changed rotation matrices in accordance with
Cu             the definition of real harmonics used in the rest of
Cu             the code (Hund's rules satisfied as indicated by orb. moment)
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   30 Apr 05 Lambrecht first created
C----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,lldau(nbas),mode,lmaxu,nsp
      double complex a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,m,l,idu(4),lmxa,is,i,j,k,ll,isp,iblu
      double complex rot(2,2)
      double complex b(2,2),c(2,2),add
      double precision s2

      s2 = 1/dsqrt(2d0)
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) /= 0) then
          is = s_site(ib)%spec
          lmxa = s_spec(is)%lmxa
          idu = s_spec(is)%idu
          do  l = 0, min(lmxa,3)
            if (idu(l+1) /= 0) then
              iblu = iblu+1
              do  isp = 1, 2
                do  m = 1, l
                  b(1,1) = a(m,m,isp,iblu)
                  b(1,2) = a(m,-m,isp,iblu)
                  b(2,1) = a(-m,m,isp,iblu)
                  b(2,2) = a(-m,-m,isp,iblu)
                  if (mode == -1) then
C rotation from spherical  to cubic basis
                    rot(1,1) = dcmplx(s2,0d0)
                    rot(1,2) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,1) = dcmplx(0d0,-s2)
C                   rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                    rot(2,1) = dcmplx(0d0,s2)
                    rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                  elseif (mode == 1) then
C rotation from cubic  to spherical  basis
                    rot(1,1) = dcmplx(s2,0d0)
C                   rot(1,2) = dcmplx(0d0,s2)
                    rot(1,2) = dcmplx(0d0,-s2)
                    rot(2,1) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                    rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                  else
                    call rx('ROTYCS: mode must be 1 or -1')
                  endif
C calculate matrix product c=rot*b*rot^+
                  do  i = 1, 2
                    do  j = 1, 2
                      add = dcmplx(0d0,0d0)
                      do  k = 1, 2
                        do  ll = 1, 2
                          add = add + rot(i,k)*b(k,ll)*dconjg(rot(j,ll))
                        enddo
                      enddo
                      c(i,j) = add
                    enddo
                  enddo
C place c in appropriate place
                  a(m,m,isp,iblu) = c(1,1)
                  a(m,-m,isp,iblu) = c(1,2)
                  a(-m,m,isp,iblu) = c(2,1)
                  a(-m,-m,isp,iblu) = c(2,2)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      end
