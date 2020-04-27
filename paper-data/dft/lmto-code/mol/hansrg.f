      subroutine hansrg(rsm,lmx,nxi,lxi,exi,r2,nrx,nr,tol,idx,wk,job,xi)
C- Vectors of smoothed Hankel functions and Gaussians
C  For each e=exi(ie), e<0, e=0 and e>0 imply respectively:
C  hankel h(rsm,exi), gaussian(rsm,exi=0), gaussian(rsm=e;exi=0)
C  exi MUST be ordered so that all hankels (e<0) come first.
C  xi is the radial part/r**l, so the solid function is
C  hl(ilm) = xi(l)*cy(ilm)*yl(ilm), so xi is radial part / r**l
C  job=1 makes true radial part.

      implicit none
      integer nrx,nr,lmx,idx(nrx,2),nxi,lxi(1),job
      double precision r2(nrx),e,xi(nrx,0:lmx,1),wk(nrx,6),
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
        if (exi(ie) < -1d-6) then
          ncut = ie
        else
          ng = ng+1
        endif
   10 continue
      if (ncut+ng /= nxi) call rx('hansrg:  energies badly ordered')

C --- need larger work array for vectorised hansr. Simplest to do this:
      call defrr(owk, (4+lmx)*nrx)
C --- Hankels and Gaussians ---
      if (ncut > 0) then
        call hansr(rsm,0,lmx,ncut,lxi,exi,r2,nrx,nr,idx,w(owk),job,xi)
      endif
      if (ng > 0) then
        do  15  ir = 1, nr
   15   wk(ir,1) = dsqrt(r2(ir))
        do  20  ig = 1, ng
          ie = ncut+ig
          rsmg = exi(ie)
          if (dabs(rsmg) < 1d-6) rsmg = rsm
          call gausr(rsmg,lxi(ie),nrx,nr,wk,xi(1,0,ie),job)
   20   continue
      endif
      call rlse(owk)

      end

C      subroutine fmain
C      implicit none
C      integer nrx,nsize,nxi,lmax
C      parameter (nrx=10000,nsize=250000,nxi=3,lmax=8)
C      double precision rofi(nrx),rsq(nrx),e,xi(nrx,0:lmax,nxi),
C     .  xi0(0:lmax),akap,
C     .  rsm,err,x,r,tol,rmsi(0:lmax),xl(0:lmax),exi(nxi)
C      integer ir,nr,lmx,l,oidx,owk,k,lxi(nxi),ie
C      real w(nsize)
C      common /static/ rofi,rsq,xi
C      common /w/ w
C      err(x,k) = (x/(xi0(k)+1d-70)-1)*1d10
C
C      call wkinit(nsize)
C      call pshprt(100)
C
C      tol = 1d-10
C      lxi(1) = 6
C      lxi(2) = 6
C      lxi(3) = 6
C      exi(1) = -1
C      exi(2) = -3
C      exi(3) = -9
C      rsm = .5d0
C      print *, 'nxi is 3'
C      print *, 'rsm,exi,tol,lxi='
C      read(*,*) rsm,exi,tol,lxi
C      call makr(rsm,nr,rofi,rsq)
C
C      akap = dsqrt(-min(exi(1),exi(2),exi(3)))
C      rsq(1) = (max(1d0,.9d0-dlog(1d10*tol)/15+akap*rsm/10)*rsm)**2
C      rsq(1) = rsq(1) - 1d-6
C      rsq(2) = rsq(1) + 2d-6
C
C      call defrr(oidx,  nr*2)
C      call defrr(owk,   nr*6)
C      call hansrg(rsm,lmax,nxi,lxi,exi,rsq,nrx,nr,tol,w(oidx),w(owk),
C     .  000,xi)
C      call rlse(oidx)
C
C      do  90, ie = 1, nxi
C      lmx = lxi(ie)
C      e = exi(ie)
C      print *, ' -------------- errors for ie=',ie,' --------------'
C      do  90  ir = 1, nr
C        r = dsqrt(rsq(ir))
C        call hnsmrx(r,e,1/(rsm+1d-25),xi0,lmx)
Cc       print 220,r,(xi0(l),l=0,lmx)
C        print 220,r,(xi(ir,l,ie),l=0,lmx),(err(xi(ir,l,ie),l),l=0,lmx)
C  220   format(8f10.6:/10x,7f10.6)
C        print *, '----'
C   90 continue
C
C      print *, 'ie for rmsi?'
C      read(*,*) ie
C      call makr2(rsm,nrx,rsq)
C      call defrr(oidx,  nrx*2)
C      call defrr(owk,   nrx*6)
C      call tm(' start')
C      call hansrg(rsm,lmax,nxi,lxi,exi,rsq,nrx,nrx,tol,w(oidx),w(owk),
C     .  000,xi)
C      call tm(' done')
C      call rlse(oidx)
C
C      lmx = lxi(ie)
C      e    = exi(ie)
C      call dpzero(rmsi,lmx+1)
C      do  190  ir = 1, nrx
C        r = dsqrt(rsq(ir))
C        call hnsmrx(r,e,1/(rsm+1d-25),xi0,lmx)
CC        if (dabs(r-3) < .001d0)
CC     .    print 220, r, (xi(ir,l),l=0,lmx),(err(xi(ir,l),l),l=0,lmx)
C        do  192 l = 0, lmx
C          if (rmsi(l) < dabs(err(xi(ir,l,ie),l))) then
C            rmsi(l) = dabs(err(xi(ir,l,ie),l))
C            xl(l) = r
C          endif
C  192   continue
C  190 continue
C      print *, 'rsmi, for energy', ie
C      print 220, (xl(l), l=0,lmx)
C      print 220, (rmsi(l), l=0,lmx)
C      call tm('exit')
C      end
C      subroutine makr2(rsm,nr,rsq)
C      implicit none
C      integer nr,ir
C      double precision rsm,rsq(1),rs
C      real ran1
C      call ran1in(1)
C
C      rs = rsm
C      if (rs == 0d0) rs = .5d0
C      do  10  ir = 1, nr
C        rsq(ir) = (10*rs*ran1())**2
C   10 continue
C      end
C      subroutine makr(rsm,nr,rofi,rsq)
C      implicit none
C      integer nr,i,ir
C      double precision rs,rsm,rsq(1),rofi(1)
C      rs = rsm
C      if (rsm < 1d-9) rs = .5d0
C      nr = 0
C      do  10  i = 0, 10
C        ir = i+1
C        rofi(ir) = dble(4*i)/10*rs
C        if (i == 0) rofi(ir) = 1.2d0
C        rsq(ir) = rofi(ir)**2
C        nr = nr+1
C   10 continue
C
C      nr = nr+1
C      rofi(nr) = 8*rs
C      rsq(nr) = rofi(nr)**2
C      nr = nr+1
C      rofi(nr) = 12*rs
C      rsq(nr) = rofi(nr)**2
C
C      nr = 3
C      rsq(3) = 13.33357d0**2
C
C      end
C
C      subroutine hnsmrx(r,e,asm,xi,lmax)
CC- Returns smoothed hankel for e<0, gaussian for e>=0
C      implicit none
C      integer lmax
C      double precision r,e,asm,xi(0:10)
C
C      if (e < -1d-6) call hansmr(r,e,asm,xi,lmax)
C      if (dabs(e) < 1d-6) call ropgau(1/asm,lmax,1,r,xi,000)
C      if (e > 1d-6) call ropgau(e,lmax,1,r,xi,000)
C      end
