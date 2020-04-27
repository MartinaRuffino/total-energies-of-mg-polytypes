      subroutine sll0sr(rsm,lmax,nlmx,nxi,lxi,exi,x,z,nrx,nr,tol,cx,
     .  rmt,llink,el,add,idx,wk,yl,isw,xi)
C- Vector of linked solid, smoothed Hankel functions (y=0)
C  isw:  1's digit: 0, interpret e >= 0 as Gaussians (hansrg)
C                   1, interpret e >= as smoothed Hankels
C       10's digit: 0, index to (lm) is (l*(l+1))/2+m+1
C                   1, index to (lm) is ((2*lmax-m+1)*m)/2+1+l
C  Note: as now implemented, assumes KE that of unsmoothed hankel at rsm
      implicit none
      integer nrx,nr,nlmx,idx(nrx,2),nxi,lxi(1),isw,llink
      double precision xi(nrx,nlmx,nxi+1),wk(nrx),exi(1),tol,rsm,cx(1),
     .  x(1),z(1),yl(nrx,1),ex(10),el,rmt,add(25,1)
      double precision adx

      integer ir,lmx,ilm,ie,l,m,lm,lmax,nlml,lmxi,iperm,ismh,lmxx,
     .  lx(10),nx

      ismh  = mod(isw,10)
      iperm = mod(isw/10,10)

C --- Radial part of smoothed hankels ---
      do  10  ir = 1, nr
   10 wk(ir) = x(ir)**2 + z(ir)**2
      lmxx = -1
      do  12  ie = 1, nxi
      lx(ie) = lxi(ie)
      ex(ie) = exi(ie)
      lmxx = max(lmxx,lxi(ie))
   12 if (lmax < lxi(ie)) call rx('sol0sr: lxi gt lmax')
      nlml = ((lmax+1)*(lmax+2))/2
      nx = nxi
      ex(nx+1) = el
      lx(nx+1) = lmxx
      if (llink == 1) nx = nx+1
      if (ismh == 0)
     .  call hansrg(rsm,nlmx-1,nx,lx,ex,wk,nrx,nr,tol,idx,yl,0,xi)
      if (ismh == 1)
     .  call hansrf(rsm,nlmx-1,nxi,lx,ex,wk,nrx,nr,tol,idx,yl,0,xi)

C --- Make ylm ---
      call ropy0m(nr,x,z,lmax,nrx,yl,wk)
      do  15  ilm = 1, nlml
      do  15  ir = 1, nr
   15 yl(ir,ilm) = cx(ilm)*yl(ir,ilm)

C --- Make yl * xi (not linked basis) ---
      if (nx > nxi) goto 100
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
      return

C --- Make yl * xi (linked basis) ---
  100 continue
      do  120  ie = 1, nxi
        lmx = lxi(ie)
        do  130  l = lmx, 0, -1
        do  130  m = l, 0, -1
          ilm = (l+1)**2-l+m
          adx = add(ilm,ie)
          lm = (l*(l+1))/2+m+1
          lmxi = lm
          if (iperm /= 0) lmxi = ((2*lmax-m+1)*m)/2 + 1 + l
          do  140  ir = 1, nr
  140     xi(ir,lmxi,ie) = yl(ir,lm)*(xi(ir,l+1,ie) + adx*xi(ir,l+1,nx))
  130   continue
  120 continue

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
C
CC      x(1) = .4d0
CC      z(1) = dsqrt(0.67d0**2 - x(1)**2)
CC      x(2) = .2d0
CC      z(2) = dsqrt(0.6904045377488563d0**2 - x(1)**2)
C      end
C      subroutine fmain
C      implicit none
C      integer nrx,nsize,nxi,lmax,nlmx
C      parameter (nrx=10000,nsize=600000,nxi=3,lmax=4,
C     .  nlmx=((lmax+1)*(lmax+2))/2)
C      double precision rofi(nrx),rsq(nrx),e,xi(nrx,nlmx,nxi+1),
C     .  xi0(0:lmax),akap,x(nrx),y(nrx),z(nrx),add(25,nxi),rsm,err,v,r,
C     .  tol,rmsi(0:lmax),xl(0:lmax),exi(nxi+1),rmt,ermx(0:10)
C      double precision cy(300),cx(300),hl(50),hx(0:10,0:10),xyz(3),el,
C     .  tl(nxi)
C      integer ir,nr,lmx,l,oidx,owk,k,lxi(nxi+1),ie,oyl,ilm,nlm,m,isw,
C     .  nhl,nxix,llink,i
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
C      if (nxi > 1) lxi(2) = lmax
C      if (nxi > 2) lxi(3) = lmax
C      exi(1) = -1
C      if (nxi > 1) exi(2) = -3
C      if (nxi > 2) exi(3) = -9
C      rsm = .5d0
C      rmt = .9
C      el = 0
C      do 77 i = 1, nxi
C   77 tl(i) = 0
C      print *, 'nxi is', nxi
C      print *, 'el,rsm,exi,lxi,rmt='
C      read(*,*) el,rsm,(exi(i),i=1,nxi),(lxi(i),i=1,nxi),rmt
C      llink = 0
C      if (el < 0) llink = 1
C      if (llink == 1) then
C        print *, 'tlink='
C        read(*,*) tl
C      endif
C      call makr(rsm,nr,x,y,z)
C      call msadd(nxi,rsm,rmt,llink,el,tl,lxi,exi,add)
C
Cc     call wkprnt(1)
C      call defrr(oidx,  nrx*2)
C      call defrr(owk,   nrx)
C      call defrr(oyl,   nrx*max(6,nlmx))
C      call wkprnt(0)
C      isw = 0
C      call sll0sr(rsm,lmax,nlmx,nxi,lxi,exi,x,z,nrx,nr,tol,cx,
C     .  rmt,llink,el,add,w(oidx),w(owk),w(oyl),isw,xi)
C      call rlse(oidx)
C
C      do  90, ie = 1, nxi
C      lmx = lxi(ie)
C      nlm = (lmx+1)**2
C      e = exi(ie)
C      print *, ' --------- errors*1e10 for ie=',ie,' --------------'
C      do  90  ir = 1, nr
C        call dpzero(ermx,10)
C        xyz(1) = x(ir)
C        xyz(2) = y(ir)
C        xyz(3) = z(ir)
C        call sllxsm(x(ir),z(ir),1/rsm,lmx,1,exi(ie),hx,rmt,el,
C     .    tl(ie),10,cx)
C        do  92  l = 0, lmx
Cc       print *, xyz
C        print 220, (xi(ir,(l*(l+1))/2+m+1,ie),m=0,l)
Cc       print 220, (hx(m,l),m=0,l)
C        print 220, (err(xi(ir,(l*(l+1))/2+m+1,ie)),m=0,l)
C        do  94  m = 0, l
C   94   ermx(m) = max(ermx(m), err(xi(ir,(l*(l+1))/2+m+1,ie)))
C   92   continue
C  220   format(8f10.6:/10x,7f10.6)
C        r = dsqrt(x(ir)**2 + z(ir)**2)
C        print *, 'done x,z,r=', sngl(x(ir)), sngl(z(ir)), sngl(r),
C     .    ' maximum errors:'
C        print 220, (ermx(l),l=0,lmx)
C        print *, '---------------'
C   90 continue
C
C      end
C      subroutine sllxsm(rho,z,asm,lmax,ne,ei,xi,rmt,el,tl,n0,cx)
CC- smoothed hankel functions, linked basis
C      implicit real*8 (a-h,p-z), integer (o)
C      dimension ei(1),xi(0:n0,0:n0,1),cx(1),
C     .   xl(200),f(0:20),fl(0:20),p(0:20),pl(0:20),tl(1)
C      call sxlm(rho,z,xl,lmax,r2)
C      r1=dsqrt(r2)
C      do 200 ie=1,ne
C      if (ei(ie) < -1d-6) call hansmr(r1,ei(ie),asm,f,lmax)
C      if (dabs(ei(ie)) < 1d-6) call ropgau(1/asm,lmax,1,r1,f,000)
C      if (ei(ie) > 1d-6) call ropgau(ei(ie),lmax,1,r1,f,000)
C      if (el < -1d-6) then
C        if (ei(ie) > -1d-6) call rx('trying to link gaussians')
C        call hansmr(r1,el,asm,fl,lmax)
C        call hansmr(rmt,ei(ie),asm,p,lmax)
C        call hansmr(rmt,el,asm,pl,lmax)
C      endif
C
C      lm=0
C      if (el > -1d-6) then
C        do 20 l=0,lmax
C        do 20 m=0,l
C        lm=lm+1
C   20   xi(m,l,ie)=f(l)*cx(lm)*xl(lm)
C      else
C        do 120 l=0,lmax
C        add = -p(l)/pl(l)*(tl(ie)-ei(ie))/(tl(ie)-el)
C        do 120 m=0,l
C        lm=lm+1
C  120   xi(m,l,ie)=(f(l) + add*fl(l))*cx(lm)*xl(lm)
C      endif
C
C  200 continue
C      end
C      subroutine msadd(nel,rsm,rmt,llink,elink,tlink,lphi,el,ceadd)
CC- Setup for linked basis
CC  ceadd: coefficients of linking energy to be added for linked.
C      implicit none
C      integer n0,nel,lphi(nel)
C      double precision elink,tlink(nel),el(nel),p(0:10),pl(0:10),
C     .  ceadd(25,5),add,rmt,rsm,wk(6),tol
C      integer lmax,ie,is,ib,l,m,ilm,ipr,idx(2),llink
C
C      call getpr(ipr)
C      ipr=51
C
CC --- Put elink to el(nel+1), (elink>=0 means no linking) ---
C      el(nel+1) = elink
C      lmax = -1
C      if (llink /= 0) then
C        do  19  ie = 1, nel
C   19   lmax = max0(lmax,lphi(ie))
C      endif
C   18 lphi(nel+1) = lmax
C
C      if (llink == 0) return
C
CC --- Stop with error if elink equals any of the elmto ---
C      do  40  ie = 1, nel
C   40 if (dabs(elink-el(ie)) < 1d-5)
C     .    call rx('elink too close to basis energy')
C
CC --- Coefficients add for linked basis ---
C      if (ipr > 50) print 346
C  346 format('  is  ie   l     ceadd')
C      tol = 1d-10
C      call hansrg(rsm,10,1,lmax,elink,rmt**2,1,1,tol,idx,wk,0,pl)
C      do  30  ie = 1, nel
C        call hansrg(rsm,10,1,lphi(ie),el(ie),rmt**2,1,1,tol,idx,wk,
C     .    0,p)
C        lmax = lphi(ie)
C        ilm = 0
C        do  34  l = 0, lmax
C        add = -p(l)/pl(l)*(tlink(ie)-el(ie))/(tlink(ie)-elink)
C        if (ipr > 50) print 345, ie,l,add
C  345   format(2i4,f12.7)
C        do  34  m = -l, l
C        ilm = ilm+1
C   34   ceadd(ilm,ie) = add
C   30 continue
C
C      end
