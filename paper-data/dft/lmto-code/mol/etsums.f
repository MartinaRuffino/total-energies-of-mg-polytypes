      subroutine etsums(rep0,rmu0,q0,repsmo,rmusmo,qsmo,rphi,
     .   reptru,rmutru,rphtru,qtru,rhoep,rhomu,rhoph,q)
c  sum the different total energy terms, printout
      implicit real*8 (a-h,p-z)
      repi=rep0-repsmo
      rmui=rmu0-rmusmo
      qi=q0-qsmo
      rhoep=repi+reptru
      rhomu=rmui+rmutru
      rhoph=rphi+rphtru
      q=qi+qtru
      if(iprint() >= 20) then
      write(6,221)
      write(6,222) 'smooth all space ', rep0,rmu0,q0
      write(6,222) 'smooth in spheres', repsmo,rmusmo,qsmo
      write(6,220) 'interstitial     ', repi,rmui,rphi,qi
      write(6,220) 'true in spheres  ', reptru,rmutru,rphtru,qtru
      write(6,220) 'total:           ', rhoep,rhomu,rhoph,q
      endif
  220 format(1x,a17,4f15.6)
  222 format(1x,a17,2f15.6,9x,' --   ',f15.6)
  221 format(/' Charge density integrals:'
     .   /26x,'rhoeps',9x,'rhomu',10x,'rhoves',10x,'q')
      end
