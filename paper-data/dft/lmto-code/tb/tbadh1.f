      subroutine tbadh1(nl,nbas,nspu,nsites,mmix,npr,delta,h0,hrs)
C- Add mixed increments to H to real-space H
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nl,nbas,nspu,nsites,mmix,npr,delta
Co Outputs:
Co   hrs: real space hamiltonian matrix
Cr Remarks
Cr   Use tbadh1 to increment H after potential mixing in dhmix
Cr   In that case delta is the diagonal
Cr   increments returned after mixing, i.e., these are the work arrays
Cr   holding all iterations of these back to mmix iterations.
C ----------------------------------------------------------------------
      implicit none
c     integer niax
c     parameter (niax=10)
C Passed Parameters
      integer nl,nbas,nspu,mmix,nsites,npr(0:1,nbas)
c    .        iax(0:niax-1,nsites)
      double precision delta(nl**2,nl**2,nbas,nspu,0:mmix+2,2),
     .                 h0(nl**2,nl**2,nsites,nspu),
     .                 hrs(nl**2,nl**2,nsites,nspu)

C Local Variables
      integer ispu,ib,j,ilm,ilmp,iprint

      if (iprint() > 40) print 100
      do  ib = 1, nbas
        do ispu = 1, nspu
          if (iprint() > 40 .and. nspu == 2) write (*,150) ispu
          j = npr(1,ib) + 1
          do  ilm = 1, nl**2
            do  ilmp = 1, nl**2

              if (ilm == ilmp .and. iprint() > 40) then
                write (*,200)
     .          ib,ilm,hrs(ilm,ilmp,j,ispu),
     .          delta(ilm,ilmp,ib,ispu,0,1),
     .          h0(ilm,ilmp,j,ispu)+delta(ilm,ilmp,ib,ispu,0,1)
              endif

              hrs(ilm,ilmp,j,ispu) = h0(ilm,ilmp,j,ispu)
     .                             + delta(ilm,ilmp,ib,ispu,0,1)
            enddo
          enddo
        enddo
      enddo

  100 format(' TBADH1: site ilm   H_in   Delta      H_out')
  150 format ('        spin ',i1)
  200 format (7x,2i4,3f10.4)
      end
