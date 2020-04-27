      subroutine pgsoscl(s_ctrl,s_pot)
C- Set SOC scaling to zero in the padded classes
      use structures
      implicit none
C ... Passed parameters
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
C ... Local parameters
      integer ib,ic

      call info2(20,1,0,' PGFASA: set SO perturbation to 0 for padded classes',0,0)
C     Set SO scale to zero for any padded classes
      do  ib  = s_ctrl%nbas+1, s_ctrl%nbasp
        ic = s_ctrl%ipc(ib)
        s_pot%socscl(ic) = 0
      enddo

      end

