      subroutine hansrf(rsm,lmx,nxi,lxi,exi,r2,nrx,nr,tol,idx,wk,job,xi)
C- Vectors of smoothed Hankel functions
C  For each e=exi(ie), e<=0, e>0 imply respectively:
C  hankel h(rsm,exi), hankel(rsm=e;exi=0)
C  exi MUST be ordered so that all hankels (e<=0) come first.
C  xi is the radial part/r**l, so the solid function is
C  hl(ilm) = xi(l)*cy(ilm)*yl(ilm), so xi is radial part / r**l
C  job=1 makes true radial part.

      implicit none
      integer nrx,nr,lmx,idx(nrx,2),nxi,lxi(1),job
      double precision r2(nrx),e,xi(nrx,0:lmx,1),wk(nrx,5),
     .  exi(1),tol,rsm,rsmg
      integer ie,ncut,ng,ig,ir,n

C ... Heap
      real w(1)
      common /w/ w
      integer owk

C --- Determine cutoff separating smoothed Hankels from Gaussians ---
      ncut = 0
      ng   = 0
      do  10  ie = 1, nxi
        if (exi(ie) < 1d-6) then
          ncut = ie
        else
          ng = ng+1
        endif
   10 continue
      if (ncut+ng /= nxi) call rx('hansrg:  energies badly ordered')

C --- need larger work array for vectorised hansr. Simplest to do this:
      call defrr(owk, (4+lmx)*nrx)
C --- Hankels of spec'd e ---
      if (ncut > 0)
     .  call hansr(rsm,0,lmx,ncut,lxi,exi,r2,nrx,nr,idx,w(owk),job,xi)

C --- Hankels of spec'd rsm ---
      if (ng > 0) then
        do  20  ig = 1, ng
          ie = ncut+ig
          rsmg = exi(ie)
          if (dabs(rsmg) < 1d-6) rsmg = rsm
          call hansr(rsmg,0,lmx,1,lxi(ie),0d0,r2,nrx,nr,idx,w(owk),job,
     .      xi(1,0,ie))
   20   continue
      endif
      call rlse(owk)

      end
