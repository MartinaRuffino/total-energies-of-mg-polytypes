      subroutine vorbydl(mode,s_site,s_spec,nbas,nl,nsp,lmaxu,lldau,pp,
     .  vorb)
C- Scale vorb by delta or 1/delta, depending on mode
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa idu
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :0 do nothing
Ci         1 scale by delta
Ci         :2 scale by 1/delta
Ci         :-1 scale by -delta
Ci         :-2 scale by -1/delta
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmaxu :dimensioning parameter for U matrix
Ci   pp    :potential parameters (atomsr.f)
Cio Inputs/Outputs
Co   vorb  :scaled by delta or 1/delta
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   17 Nov 07 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nbas,nl,nsp,lldau(nbas),lmaxu
      double precision pp(6,nl,nsp,*)
      integer lmxa,idu(4)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer ib,is,ic,l,m1,m2
      integer iblu
      double precision xx,enu,calp,gam,alp,dela,delg,fac

      if (mode == 0) return
      call sanrg(.true.,mode,-2,2,'vorbydl','mode')
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) /= 0) then
          is = s_site(ib)%spec
          ic = s_site(ib)%class
          lmxa = s_spec(is)%lmxa
          idu = s_spec(is)%idu
          do  l = 0, min(lmxa,3)
          if (idu(l+1) /= 0) then
            iblu = iblu+1

            do  is = 1, nsp
              enu = pp(1,l+1,is,ic)
              calp = pp(2,l+1,is,ic)
              gam = pp(5,l+1,is,ic)
              alp = pp(6,l+1,is,ic)
              dela = pp(3,l+1,is,ic)**2
              xx = 1 + (calp-enu)*(gam-alp)/dela
              delg = dela*xx**2
              fac = delg
              if (iabs(mode) == 2) fac = 1/delg
              if (mode < 0) fac = -fac
C             fac = 1
              do  m2 = -l, l
              do  m1 = -l, l
                vorb(m1,m2,is,iblu) = vorb(m1,m2,is,iblu)*fac
              enddo
              enddo
            enddo

          endif
          enddo
        endif
      enddo

C      call zprm('vorb after scaling',2,vorb,
C     .  (2*lmaxu+1)**2,(2*lmaxu+1)**2,nsp*iblu)


      end
