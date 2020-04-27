      subroutine pgbasp(nbas,npadl,npadr,bas,plat,platl,platr,basp)
C- Pads basis to accommodate extra bulk layers on LHS and RHS.
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   npadl
Ci   npadr
Ci   bas   :basis vectors, in units of alat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   platl
Ci   platr
Co Outputs
Co   basp  :padded basis
C ----------------------------------------------------------------
      implicit none
      integer nbas,npadl,npadr
      double precision plat(3,3),bas(3,nbas),basp(3,*),
     .  platl(3,3),platr(3,3)
C ... Local parameters
      double precision wk(3),pgfn(3),xx,ddot
      integer ib,nbasp

C     call prmx('original basis',bas,3,3,nbas)

      call dpcopy(bas,basp,1,3*nbas,1d0)

C --- make pgfn ---
      call cross(plat,plat(1,2),pgfn)
      xx = ddot(3,pgfn,1,plat(1,3),1)
      if (xx < 0) call dscal(3,-1d0,pgfn,1)

C --- Padding on LHS ---
      nbasp = nbas
      xx = ddot(3,pgfn,1,platl(1,3),1)
      call dpcopy(platl(1,3),wk,1,3,dsign(1d0,xx))
      do  ib = 1, npadl
        call dpcopy(bas(1,ib),basp(1,nbasp+1),1,3,1d0)
        call dpadd(basp(1,nbasp+1),wk,1,3,-1d0)
        nbasp = nbasp+1
      enddo

C --- Padding on RHS ---
      xx = ddot(3,pgfn,1,platr(1,3),1)
      call dpcopy(platr(1,3),wk,1,3,dsign(1d0,xx))
      do  ib = nbas-npadr+1, nbas
        call dpcopy(bas(1,ib),basp(1,nbasp+1),1,3,1d0)
        call dpadd(basp(1,nbasp+1),wk,1,3,1d0)
        nbasp = nbasp+1
      enddo

C      call prmx('padded basis',basp,3,3,nbasp)
C      call fopna('out',66,0)
C      call dfdump(basp,3*nbasp,-66)
C      stop
C      do  50  ib = 1, nbasp
C        print 333, ib, basp(1,ib), basp(2,ib), basp(3,ib)
C  333   format(i4,3f12.6)
C   50 continue

      end
