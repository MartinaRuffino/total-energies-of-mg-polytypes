      subroutine sol0sr(rsm,lmax,nlmx,nxi,lxi,exi,x,z,nrx,nr,tol,cx,
     .  idx,wk,yl,isw,xi)
C- Vector of solid, smoothed Hankel functions, set of negative e's (y=0)
C  isw:  1's digit: 0, interpret e >= 0 as Gaussians (hansrg)
C                   1, interpret e >= as smoothed Hankels
C       10's digit: 0, index to (lm) is (l*(l+1))/2+m+1
C                   1, index to (lm) is ((2*lmax-m+1)*m)/2+1+l
      implicit none
      integer nrx,nr,nlmx,idx(nrx,2),nxi,lxi(1),isw
      double precision xi(nrx,nlmx,1),wk(nrx),exi(1),tol,rsm,cx(1),
     .  x(1),z(1),yl(nrx,1)
      integer ir,ll,lmx,ilm,ie,l,m,lm,lmax,nlml,lmxi,iperm,ismh

      ismh  = mod(isw,10)
      iperm = mod(isw/10,10)

C --- Make r**2 ---
      do  10  ir = 1, nr
   10 wk(ir) = x(ir)**2 + z(ir)**2
      do  12  ie = 1, nxi
   12 if (lmax < lxi(ie)) call rx('sol0sr: lxi gt lmax')
      nlml = ((lmax+1)*(lmax+2))/2

      if (ismh == 0)
     .  call hansrg(rsm,nlmx-1,nxi,lxi,exi,wk,nrx,nr,tol,idx,yl,0,xi)
      if (ismh == 1)
     .  call hansrf(rsm,nlmx-1,nxi,lxi,exi,wk,nrx,nr,tol,idx,yl,0,xi)
      call ropy0m(nr,x,z,lmax,nrx,yl,wk)

      do  15  ilm = 1, nlml
      do  15  ir = 1, nr
   15 yl(ir,ilm) = cx(ilm)*yl(ir,ilm)

C --- Make yl * xi ---
      do  20  ie = 1, nxi
        lmx = lxi(ie)
        do  30  l = lmx, 0, -1
        do  30  m = l, 0, -1
          lm = (l*(l+1))/2+m+1
          lmxi = lm
          if (iperm /= 0) lmxi = ((2*lmax-m+1)*m)/2 + 1 + l
          do  40  ir = 1, nr
   40     xi(ir,lmxi,ie) = yl(ir,lm)*xi(ir,l+1,ie)
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
C        x(i) = abs((ran1()-.5d0)*5*rs)
C        y(i) = (ran1()-.5d0)*5*rs*0
C        z(i) = (ran1()-.5d0)*5*rs
C   10 continue
C      end
C      subroutine fmain
C      implicit none
C      integer nrx,nsize,nxi,lmax,nlmx
C      parameter (nrx=10000,nsize=600000,nxi=3,lmax=4,
C     .  nlmx=((lmax+1)*(lmax+2))/2)
C      double precision rofi(nrx),rsq(nrx),e,xi(nrx,nlmx,nxi),
C     .  xi0(0:lmax),akap,x(nrx),y(nrx),z(nrx),
C     .  rsm,err,v,r,tol,rmsi(0:lmax),xl(0:lmax),exi(nxi)
C      double precision cy(300),cx(300),hl(50),hx(0:10,0:10),xyz(3)
C      integer ir,nr,lmx,l,oidx,owk,k,lxi(nxi),ie,oyl,ilm,nlm,m,isw
C      real w(nsize)
C      common /static/ rofi,rsq,xi
C      common /w/ w
C      err(v) = (hx(m,l)/v-1)*1d10
C
C      call wkinit(nsize)
C      call pshprt(100)
C      call sylmnc(cy,12)
C      call sxlmnc(cx,10)
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
C      isw = 0
C      call sol0sr(rsm,lmax,nlmx,nxi,lxi,exi,x,z,nrx,nr,tol,cx,
C     .  w(oidx),w(owk),w(oyl),isw,xi)
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
C        call solxsm(x(ir),z(ir),1/rsm,lmx,1,exi(ie),hx,10,cx)
C        do  92  l = 0, lmx
Cc       print *, xyz
C        print 220, (xi(ir,(l*(l+1))/2+m+1,ie),m=0,l)
Cc       print 220, (hx(m,l),m=0,l)
C        print 220, (err(xi(ir,(l*(l+1))/2+m+1,ie)),m=0,l)
C   92   continue
C  220   format(8f10.6:/10x,7f10.6)
C        print *, '----'
C   90 continue
C
C      end
C      SUBROUTINE SOLXSM(RHO,Z,ASM,LMAX,NE,EI,XI,N0,CX)
CC  MAKES SMOOTHED HANKEL FUNCTIONS
C      IMPLICIT REAL*8 (A-H,P-Z), INTEGER (o)
C      DIMENSION EI(1),XI(0:N0,0:N0,1),CX(1),
C     .   XL(200),F(0:20)
C      CALL SXLM(RHO,Z,XL,LMAX,R2)
C      R1=DSQRT(R2)
C      DO 20 IE=1,NE
C      if (ei(ie) < -1d-6) call hansmr(r1,ei(ie),asm,f,lmax)
C      if (dabs(ei(ie)) < 1d-6) call ropgau(1/asm,lmax,1,r1,f,000)
C      if (ei(ie) > 1d-6) call ropgau(ei(ie),lmax,1,r1,f,000)
C
CC      if (ei(ie) > -1d-6) then
CC        print 345, (f(i), i=0,lmax)
CC  345   format(6f12.6)
CC        pause
CC      endif
C
C      LM=0
C      DO 20 L=0,LMAX
C      DO 20 M=0,L
C      LM=LM+1
C  20  XI(M,L,IE)=F(L)*CX(LM)*XL(LM)
C
C      END
