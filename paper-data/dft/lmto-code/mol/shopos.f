      subroutine shopos(nbas,nbasp,pos,mpole,dpole,atid)
C- Display table of atom positions and point multipoles
      implicit none
      integer nbas,nbasp
      double precision pos(3,nbasp),mpole(nbasp),dpole(3,nbasp)
      character*8 atid(1)

      integer ib,i,iprint

      if (iprint() < 30) return

      print *
      if (nbasp > nbas) then
        write(*,50)
      else
        write(*,51)
      endif
      do  ib = 1, nbas
        write (*,100) atid(ib), (pos(i,ib),i=1,3)
      enddo
      if (nbasp > nbas) then
        do  ib = nbas+1, nbasp
          write (*,200) atid(ib), (pos(i,ib),i=1,3), mpole(ib),
     .                        (dpole(i,ib),i=1,3)
        enddo
      endif
      print *

   50 format (' Species',17x,'position',20x,'Q',11x,
     .        'px',10x,'py',10x,'pz')
   51 format (' Species',17x,'position')
  100 format (4x,a4,1x,3f12.6)
  200 format (4x,a4,1x,7f12.6)
      end
