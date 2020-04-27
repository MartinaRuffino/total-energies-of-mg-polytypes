      subroutine spvls2(nlml,lmxa,nrad,nel,el,rad,nbas,vl,indxcg,
     .  jcg,cg,pjj,npjj,ioffp0)
C- Kinetic,potential and overlap perturbation matrices for smoothed pot
C  vl is smoothed potential, scaled by weights
      implicit none
      integer lmxa,nrad,nel,nbas,
     . lphi(1),indxcg(1),jcg(1),nlml,ioffp0,npjj
      double precision vl(nrad,1),cg(1),el(1),pjj(npjj,1),rad(1)
      integer ie,iv,je,job,lsp,nlma,npan,nsp,oy,oh,obi,obj,osum
      real w(1)
      common /w/ w

C --- Setup ---
      nsp = lsp()+1
      npan = 1
      if (nsp /= 1) call rx('spvls not set up for spin pol')
      if (npan /= 1) call rx('spvls not set up for tpan')
      nlma = (lmxa+1)**2
      call defrr(oy,    nrad)
      call defrr(oh,    nrad)
      call defrr(obi,   nrad*(lmxa+1))
      call defrr(obj,   nrad*(lmxa+1))
      call defrr(osum,  nlml*(lmxa+1)**2)

C --- Make pjj = (J V J) for each pair of energies ---
      iv = 0
      do  20  ie = 1, nel
        call ropj(rad,el(ie),lmxa,w(oy),w(oh),w(obi),nrad)
        do  22  je = ie, nel
          iv = iv+1
          call ropj(rad,el(je),lmxa,w(oy),w(oh),w(obj),nrad)
          job = 2
          if (ie == je) job = 1
          call pertmw(nlma,nlma,w(obi),w(obj),vl,nlml,job,
     .      w(osum),lmxa+1,nrad,cg,jcg,indxcg,pjj(ioffp0+1,iv))
   22   continue
   20 continue

      call rlse(oy)
      end
