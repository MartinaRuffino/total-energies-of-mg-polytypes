      subroutine solhsr(rsm,nlmx,nxi,lxi,exi,x,y,z,nrx,nr,tol,cy,
     .  idx,wk,yl,xi)
C- Vector of solid, smoothed Hankel functions, set of negative e's
      implicit none
      integer nrx,nr,nlmx,idx(nrx,2),nxi,lxi(1)
      double precision xi(nrx,nlmx,1),wk(nrx),exi(1),tol,rsm,cy(1),
     .  x(1),y(1),z(1),yl(nrx,1)
      integer ir,ll,lmx,ilm,ie,l,m

C --- Make r**2 ---
      do  10  ir = 1, nr
   10 wk(ir) = x(ir)**2 + y(ir)**2 + z(ir)**2

      call hansrg(rsm,nlmx-1,nxi,lxi,exi,wk,nrx,nr,tol,idx,yl,0,xi)
      call ropylm(nr,x,y,z,ll(nlmx),nrx,yl,wk)

      do  15  ilm = 1, nlmx
      do  15  ir = 1, nr
   15 yl(ir,ilm) = cy(ilm)*yl(ir,ilm)

C --- Make yl * xi ---
      do  20  ie = 1, nxi
        lmx = lxi(ie)
        ilm = (lmx+1)**2
        do  30  l = lmx, 0, -1
        do  30  m = -l, l
          do  40  ir = 1, nr
   40     xi(ir,ilm,ie) = yl(ir,ilm)*xi(ir,l+1,ie)
          ilm = ilm-1
   30   continue
   20 continue

      end
C      subroutine makr(rsm,nr,x,y,z)
C      implicit none
C      integer nr,i,ir
C      double precision rs,rsm,x(1),y(1),z(1)
C      real ran1
C      rs = rsm
C      if (rsm < 1d-9) rs = .5d0
C      call ran1in(1)
C      nr = 5
C      do  10  i = 1, nr
C        ir = i+1
C        x(i) = (ran1()-.5d0)*5*rs
C        y(i) = (ran1()-.5d0)*5*rs
C        z(i) = (ran1()-.5d0)*5*rs
C   10 continue
C
C      end
C      subroutine fmain
C      implicit none
C      integer nrx,nsize,nxi,lmax,nlmx
C      parameter (nrx=10000,nsize=600000,nxi=3,lmax=2,nlmx=(lmax+1)**2)
C      double precision rofi(nrx),rsq(nrx),e,xi(nrx,nlmx,nxi),
C     .  xi0(0:lmax),akap,x(nrx),y(nrx),z(nrx),
C     .  rsm,err,v,r,tol,rmsi(0:lmax),xl(0:lmax),exi(nxi)
C      double precision cy(300),hl(50),xyz(3)
C      integer ir,nr,lmx,l,oidx,owk,k,lxi(nxi),ie,oyl,ilm,nlm
C      real w(nsize)
C      common /static/ rofi,rsq,xi
C      common /w/ w
C      err(v,k) = (v/(hl(k)+1d-70)-1)*1d10
C
C      call wkinit(nsize)
C      call pshprt(100)
C      call sylmnc(cy,12)
C
C      tol = 1d-10
C      lxi(1) = lmax
C      lxi(2) = lmax
C      lxi(3) = lmax
C      exi(1) = -1
C      exi(2) = -3
C      exi(3) = -9
C      rsm = .5d0
C      print *, 'nxi is 3'
C      print *, 'rsm,exi,tol,lxi='
C      read(*,*) rsm,exi,tol,lxi
C      call makr(rsm,nr,x,y,z)
C
Cc     call wkprnt(1)
C      call defrr(oidx,  nrx*2)
C      call defrr(owk,   nrx)
C      call defrr(oyl,   nrx*max(6,nlmx))
C      call wkprnt(0)
C      call solhsr(rsm,nlmx,nxi,lxi,exi,x,y,z,nrx,nr,tol,cy,
C     .  w(oidx),w(owk),w(oyl),xi)
C      call rlse(oidx)
C
C      do  90, ie = 1, nxi
C      lmx = lxi(ie)
C      nlm = (lmx+1)**2
C      e = exi(ie)
C      print *, ' -------------- errors for ie=',ie,' --------------'
C      do  90  ir = 1, nr
C        xyz(1) = x(ir)
C        xyz(2) = y(ir)
C        xyz(3) = z(ir)
C        call solhsm(xyz,exi(ie),1/rsm,lmx,hl,cy)
C        print 220, xyz
C        print 220, (xi(ir,ilm,ie),ilm=1,nlm),
C        print 220, (err(xi(ir,ilm,ie),ilm),ilm=1,nlm)
C  220   format(8f10.6:/10x,7f10.6)
C        print *, '----'
C   90 continue
C
C      end
C      SUBROUTINE SOLHSM(R,E,A,LMAX,HL,CY)
CC  REAL SMOOTHED SOLID HANKEL FUNCTIONS, NEGATIVE ENERGIES
C      IMPLICIT REAL*8 (A-H,P-Z), INTEGER(o)
C      DIMENSION R(3),HL(1),BL(1),CY(1),PSI(0:20)
C      CALL SYLM(R,HL,LMAX,R2)
C      R1=DSQRT(R(1)**2+R(2)**2+R(3)**2)
C
C      if (e < -1d-6) call hansmr(r1,e,a,psi,lmax)
C      if (dabs(e) < 1d-6) call ropgau(1/a,lmax,1,r1,psi,000)
C      if (e > 1d-6) call ropgau(e,lmax,1,r1,psi,000)
C      ILM=0
C      DO 10 L=0,LMAX
C      NM=2*L+1
C      DO 10 M=1,NM
C      ILM=ILM+1
C  10  HL(ILM)=PSI(L)*CY(ILM)*HL(ILM)
C      RETURN
C      END
