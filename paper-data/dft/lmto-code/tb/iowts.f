      subroutine iowts(nsp,nevx,nbmax,eband,wtkb,efermi,qval)
C- write out eigenvalue weights for case MOL=T
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nsp,nevx,wtkb
Co Outputs:
Co   weights written to disc
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters
      integer nevx,nsp,nbmax
      double precision eband(nbmax,nsp),wtkb(nevx,nsp),efermi,qval(2)
C Local Variables
      integer i,isp,ifi,fopn

      ifi = fopn('EV')
      if (nsp == 1) then
        call awrit1('Evals, weights: dimension %i',' ',56,ifi,nevx)
        write (ifi,'(2f15.9)') (eband(i,1),wtkb(i,1), i = 1, nevx)
        call awrit2('Fermi energy %d, %d electrons',' ',128,ifi,efermi,
     .              qval(1))
      else
        call awrit1('Evals weights (up, down spin): dimension %i',' ',
     .              128,ifi,nevx)
        write (ifi,'(4f15.9)') (eband(i,1), wtkb(i,1), eband(i,2),
     .                          wtkb(i,2), i = 1, nevx)
       call awrit3('Fermi energy %d, %d electrons, magnetic moment: %d',
     .              ' ',256,ifi,efermi,qval(1),qval(2))
      endif
      end
