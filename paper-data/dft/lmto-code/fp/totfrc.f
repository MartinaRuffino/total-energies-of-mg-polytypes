      subroutine totfrc(nbas,s_site,s_lat,leks,fes1,fes2,dfhf,f,fmax,dfmax)
C- Add together and print contributions to force
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     force
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:symgr istab
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   leks: :<0 no forces calculated
Ci         : 0 Harris forces only
Ci         :>0 also print HKS forces
Ci         :>1 use HKS forces instead of HF forces
Ci   fes1  :contribution to HF forces from estat + xc potential
Ci   fes2  :contribution to KS forces from estat + xc potential
Ci   dfhf  :2nd order corr. HF forces from ansatz density shift (dfrce)
Cio Inputs/Outputs
Cio     f  :On input, f is the contribution to force from eigval sum
Cio        :On output, f is the total force
Co Outputs
Ci   fmax  :largest force
Ci   dfmax :largest correction term
Cr Remarks
Cu Updates
Cu   21 Sep 17 Returns c*fmax and c*dfmax
Cu   03 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 Apr 03 Prints out max correction to Harris force
Cu   30 May 00 Adapted from nfp totforce.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,leks
      double precision f(3,nbas),fes1(3,nbas),fes2(3,nbas),dfhf(3,nbas),fmax,dfmax
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer stdo,stdl,ipr,ipl,ib,m,ng,ibmx,isw
      double precision fev1(3),fev2(3),fhar(3),fks(3),c
      procedure(integer) :: lgunit
      procedure(real(8)) :: ddot

      stdo = lgunit(1)
      stdl = lgunit(2)
      if (leks < 0) return
      c = 1000d0
      call getpr(ipr)
      ipl = 1
      fmax = -1
      dfmax = -1

      if (ipr >= 30) then
        if (ddot(3*nbas,dfhf,1,dfhf,1) /= 0) then
          write(stdo,'(/'' Forces, with eigenvalue correction'')')
        else
          write(stdo,'(/''Forces:'')')
        endif
        write(stdo,201)
  201   format('  ib',11x,'estatic',18x,'eigval',20x,'total')
      endif
      do  ib = 1, nbas
        do  m = 1, 3
          fev1(m) = f(m,ib) + dfhf(m,ib)
          fev2(m) = f(m,ib)
          fhar(m) = fev1(m) + fes1(m,ib)
          f(m,ib) = fhar(m)
        enddo
        if (dsqrt(ddot(3,fhar,1,fhar,1)) > fmax) then
          ibmx = ib
          fmax = dsqrt(ddot(3,fhar,1,fhar,1))
        endif
        if (dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1)) > dfmax) then
          ibmx = ib
          dfmax = dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1))
        endif
C       if (ipr >= 30) write(stdo,200) ib,(c*fes1(m,ib),m=1,3),(c*fev1(m),m=1,3),(c*fhar(m),m=1,3)
        call info5(30,0,0,'%,4i%3;8,2D %3;8,2D %3;8,2D',ib,c*fes1(:,ib),c*fev1,c*fhar,0)
  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2)
        if (leks. eq. 0 .and. ipl > 1 .and. ipr >= 30) write(stdl,710) ib,(c*fhar(m),m=1,3)
        s_site(ib)%force(1:3) = f(1:3,ib)
  710   format('fp ib',i4,'  fh ',3f8.2,2x,3f8.2)

        if (leks >= 1) then
          do  m = 1, 3
            fks(m)  = fev2(m) + fes2(m,ib)
            if (leks > 1) f(m,ib) = fks(m)
            if (leks > 1) s_site(ib)%force = fks
          enddo
C         if (ipr > 40) write(stdo,210) (c*fes2(m,ib),m=1,3),(c*fev2(m),m=1,3),(c*fks(m),m=1,3)
          call info5(41,0,0,'  KS%3;8,2D %3;8,2D %3;8,2D',c*fes2(:,ib),c*fev2,c*fks,0,0)
  210     format('  KS',3f8.2,1x,3f8.2,1x,3f8.2)
          if (ipl > 1 .and. ipr >= 30) write (stdl,711) ib,(c*fhar(m),m=1,3),(c*fks(m),m=1,3)
  711     format('fp ib',i4,'  fh ',3f8.2,'   fks ',3f8.2)
        endif
      enddo
      call info5(10,0,0,' Maximum Harris force = %;3g mRy/au (site %i)'
     .  //'%?#n#  Max eval correction = %;1d##',c*fmax,ibmx,isw(dfmax > 0),c*dfmax,0)
      fmax = c*fmax; dfmax = c*dfmax

C     Symmetrize forces to machine precision
      ng = s_lat%nsgrp
      if (ng > 1) then
        call info(30,1,0,' Symmetrize forces ...',0,0)
        call symfor(nbas,1,s_lat%symgr,ng,s_lat%istab,f)
      endif

      end
