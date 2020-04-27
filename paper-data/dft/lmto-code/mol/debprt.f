      subroutine debprt(string,nsp,nr,rps)
C- Print an array for debugging
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character string
      integer nsp,nr
      double precision rps(nr)
      integer ip

      stop 'np not defined'

C       write(6,'(1x,3a)') ' ************ ',string,' ****'
C       if(nsp == 1) then
C        write(6,*) ' nr =',nr
C        do ip=1,np,10
C          write(6,'(i6,100g12.4)') ip,(rps(ir),ir=1,nr,10)
C        enddo
C        enddo
C       else
C        do isp = 1,nsp
C          write(6,*) ' nr, spin =',nr,isp
C          do ip=1,np,10
C            write(6,'(i6,100g12.4)') ip,(rps(ir+nr*(isp-1)),ir=1,nr,10)
C          enddo
C        enddo
C       endif
       end
