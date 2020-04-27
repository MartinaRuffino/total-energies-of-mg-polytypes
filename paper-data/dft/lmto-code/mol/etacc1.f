c --------- etacc1: init energy terms to zero --------------
      subroutine etacc1(repsmo,rmusmo,rphsmo,qsmo,reptru,rmutru,
     .   rphtru,qtru)
      implicit real*8 (a-h,p-z)
      repsmo=0d0
      rmusmo=0d0
      rphsmo=0d0
      qsmo=0d0
      reptru=0d0
      rmutru=0d0
      rphtru=0d0
      qtru=0d0
      end
c --------- etacc2: add to accumulated energies ---------
      subroutine etacc2(rep1,rmu1,rph1,q1,repsmo,rmusmo,rphsmo,qsmo)
      implicit real*8 (a-h,p-z)
      repsmo=repsmo+rep1
      rmusmo=rmusmo+rmu1
      rphsmo=rphsmo+rph1
      qsmo=qsmo+q1
      end
c --------- etacc3: make final energy terms and printout ----------
      subroutine etacc3(repsmo,rmusmo,rphsmo,qsmo,reptru,rmutru,
     .   rphtru,qtru,rhoep,rhomu,rhoph,q)
      implicit real*8 (a-h,p-z)
      write(6,221)
      write(6,220) 'smooth all space ', rhoep,rhomu,rhoph,q
      write(6,220) 'smooth in spheres', repsmo,rmusmo,rphsmo,qsmo
      write(6,220) 'true in spheres  ', reptru,rmutru,rphtru,qtru
      rhoep=rhoep-repsmo
      rhomu=rhomu-rmusmo
      rhoph=rhoph-rphsmo
      q=q-qsmo
      write(6,220) 'interstitial     ', rhoep,rhomu,rhoph,q
      rhoep=rhoep+reptru
      rhomu=rhomu+rmutru
      rhoph=rhoph+rphtru
      q=q+qtru
      write(6,220) 'total:           ', rhoep,rhomu,rhoph,q
  220 format(1x,a17,4f15.6)
  221 format(/' Charge density integrals:'
     .   /26x,'rhoeps',9x,'rhomu',10x,'rhoves',10x,'q')
      end
