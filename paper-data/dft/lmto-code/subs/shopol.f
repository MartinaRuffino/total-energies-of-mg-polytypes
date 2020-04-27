      subroutine shopol(nl,nclass,nsp,s_spec,ics,pp)
C- Find poles for third order potential functions
C-----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nl,nbas,nsp,isp,ipc,pp (in gamma representation)
Co Outputs
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nl,nsp,nclass,ics(nclass)
      double precision pp(6,nl,nsp,*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer l,ic,isp,ipr,i,lmxa,is
      double precision oalpha,pgam,enu,C,gam,del,pa,alp,xx
      double complex zz,phold(3),zmc,guess,update

      call getpr(ipr)
      ipr = 55
      if (ipr >= 55) then
        if (nsp == 1) then
          print '(''  ib   l   Poles P_l ...'')'
        else
          print '(''  ib spin  Poles P_l(ib) ...'')'
        endif
      endif

      do  isp = 1, nsp
        do  ic = 1, nclass
          is = ics(ic)
          lmxa = s_spec(is)%lmxa
          do  l = 0, nl-1
            if (l > lmxa) cycle
            enu = pp(1,l+1,isp,ic)
            C   = pp(2,l+1,isp,ic)
            del = pp(3,l+1,isp,ic)**2
            pa  = pp(4,l+1,isp,ic)
            gam = pp(5,l+1,isp,ic)
            alp = pp(6,l+1,isp,ic)
            pgam = pa - oalpha(enu,C,del,alp,gam)**2

            guess = C
C           guess = (-1.045d0, .002d0)
            i = 1
    1       continue
            zz = guess + pgam * (guess - enu)**3
            zmc = zz - C
            update = zmc/(1+pgam*3*(guess-enu)**2)
            guess = guess - update
            if (cdabs(update) > 1d-10) goto 1
            phold(i) = guess
            if (i == 1) then
              i = 2
              xx = (C-enu)**2*pgam
              guess = dcmplx(enu-(C-enu)*(1-xx)/2,
     .                       sqrt(1/pgam)*(1+3*xx/8))
              goto 1
            endif
            if (i == 2) then
              i = 3
              guess = dconjg(guess)
              goto 1
            endif
            if (ipr >= 55 .and. nsp == 1) then
              print 333,  ic, l, (phold(i), i=1,3)
            elseif (ipr >= 55) then
              print 334,  ic, isp, l, (phold(i), i=1,3)
            endif
          enddo
        enddo
      enddo
  333 format(2i4,10(1x,2f10.5,2x))
  334 format(3i4,10(1x,2f10.5,1x))
      end
