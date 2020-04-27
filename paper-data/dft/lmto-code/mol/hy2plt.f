      subroutine hy2plt(lx1,ex1,nx1,cus1,npwr1,lx2,ex2,nx2,cus2,npwr2,
     .  dr,r1,r2,rsm1,rsm2,e1,e2,c1,c2,nf1,nf2,nlm1,nlm2)
C- Compare rotated tc-fit wave function product
      implicit none
      integer lx1(1),lx2(1)
      integer nf1,nf2,nlm1,nlm2,nx1,nx2,npwr1,npwr2
      double precision rmax(2),dr(3),r1,r2,rsm1,rsm2,
     .  tau(3,2),xv(3),yv(3),zv(3),cpar(3,2),cus1(npwr1,1),cus2(npwr2,1)
      double precision ex1(1),ex2(1),e1(1),e2(1),
     .  c1(nf1,nlm1,nlm2),c2(nf2,nlm2,nlm1)
      integer i,ib,ilm1,ilm2,ix,iy,j,l,l1,l2,ll,lp1,lp2,m,
     .  nx,nx0,ny,ny0,ilm0,jlm0,ipr,oph1,oph2,oidx,oxi,owk,oyl,
     .  ip,iprint,lmax,nxi,lxi(10),nlm,nlmp1,nlmp2,np,npmx,
     .  np1,np2,ip1,ip2,jp,oul,osl,ophr1,ophr2,i1,i2
      parameter(  nx0=25, ny0=25, npmx=nx0*ny0)
      integer idx(npmx)
      real prod(nx0,ny0),fit(nx0,ny0),xszcm,yszcm
      double precision dh,h1,h2,sp3,prd,sum,x1,x2,xnrm,
     .  xx,xxx,y1,y2,ynrm,yyy,z0,zdist,bot,top,
     .  cy(100),ph1(100),ph2(100),p(3),p1(3),p2(3),rad(npmx),
     .  fac1(0:10),fac2(0:10),fac,xi1(500),xi2(500),
     .  xyz(npmx,3),rsmp1,rsmp2
      real w(1)
      common /w/ w

C --- Setup ---
      call getpr(ipr)
      call sylmnc(cy,8)
      lp1 = ll(nlm1)
      lp2 = ll(nlm2)
      do  1  m = 1, 3
      tau(m,1) = -dr(m)*0.5d0
    1 tau(m,2) =  dr(m)*0.5d0
      rmax(1) = r1
      rmax(2) = r2
      write(6,338) tau
  338 format(' tau1=',3f9.4,'   tau2=',3f9.4)
C --- Get max values on sphere to normalize phi1,phi2 ---
      call hrmsav(r1,lp1,e1,fac1)
      call hrmsav(r2,lp2,e2,fac2)
      if (ipr >= 30) print 330, 1,(fac1(l),l=0,lp1)
      if (ipr >= 30) print 330, 2,(fac2(l),l=0,lp2)
  330 format(' hy2plt: rms avg phi at r',i1,': ',6f9.4)
C --- Define plane ---
      write(6,*) 'For plane enter xv,yv,xmin,xmax,ymin,ymax,z0'
      read(5,*)  xv,yv,x1,x2,y1,y2,z0
      xnrm = dsqrt(sp3(xv,xv))
      do  10  m = 1, 3
   10 xv(m) = xv(m)/xnrm
      xx = sp3(xv,yv)
      do  11  m = 1, 3
   11 yv(m) = yv(m) - xx*xv(m)
      ynrm = dsqrt(sp3(yv,yv))
      do  12  m = 1, 3
   12 yv(m) = yv(m)/ynrm
      call cross(xv,yv,zv)
      write(6,882) xv,yv,zv
  882 format(' xv=',3f7.3,'   yv=',3f7.3,'   zv=',3f7.3)
      do  13  ib = 1, 2
      cpar(1,ib) = sp3(xv,tau(1,ib))
      cpar(2,ib) = sp3(yv,tau(1,ib))
      zdist = sp3(zv,tau(1,ib)) - z0
      cpar(3,ib) = .01
      xxx = rmax(ib)**2 - zdist**2
      if (xxx > 0d0) cpar(3,ib) = dsqrt(xxx)
   13 write(6,997) ib,cpar(1,ib),cpar(2,ib),cpar(3,ib)
  997 format(' ib=',i5,'   cpar=',3f9.5)
c --- Select pictures ---
      write(6,*) 'Enter nx,ny,  pairs ilm,jlm'
      read (5,*)  nx,ny, ilm0,jlm0
      l1 = ll(ilm0)
      l2 = ll(jlm0)
      fac = 1d0/(fac1(l1)*fac2(l2))
      bot =  1d10
      top = -1d10
c --- Generate mesh, excluding points inside sphere ---
      np = nx*ny
      np1 = 0
      np2 = 0
      np = 0
      do  50  ix = 1, nx
      do  50  iy = 1, ny
      xxx = x1
      yyy = y1
      if (nx > 1) xxx = x1 + (ix-1)*(x2-x1)/(nx-1d0)
      if (ny > 1) yyy = y1 + (iy-1)*(y2-y1)/(ny-1d0)
      do  52  m = 1, 3
      p(m) = xxx*xv(m) + yyy*yv(m) + z0*zv(m)
      p1(m) = p(m)-tau(m,1)
      p2(m) = p(m)-tau(m,2)
   52 continue
      if (p1(1)**2+p1(2)**2+p1(3)**2 < r1**2) then
        np1 = np1+1
c        print 356, '1 in np1:', ix,iy,p,yyy,y1,y2
      elseif (p2(1)**2+p2(2)**2+p2(3)**2 < r2**2) then
        np2 = np2+1
c       print 356, '1 in np2:', ix,iy,p,yyy,y1,y2
      else
        np = np+1
      endif
   50 continue
c --- Generate mesh, appending points inside spheres to top ---
      ip = 0
      ip1 = np
      ip2 = np+np1
      np = 0
      do  55  ix = 1, nx
      do  55  iy = 1, ny
      np = np+1
      xxx = x1
      yyy = y1
      if (nx > 1) xxx = x1 + (ix-1)*(x2-x1)/(nx-1d0)
      if (ny > 1) yyy = y1 + (iy-1)*(y2-y1)/(ny-1d0)
      do  56  m = 1, 3
      p(m) = xxx*xv(m) + yyy*yv(m) + z0*zv(m)
      p1(m) = p(m)-tau(m,1)
      p2(m) = p(m)-tau(m,2)
   56 continue
      if (p1(1)**2+p1(2)**2+p1(3)**2 < r1**2) then
C        print 356, '2 in np1:', ix,iy,p,yyy,y1,y2
C  356   format(a10,2i3,6f12.5)
        ip1 = ip1+1
        jp = ip1
        rad(jp) = dsqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
      elseif (p2(1)**2+p2(2)**2+p2(3)**2 < r2**2) then
C       print 356, '2 in np2:', ix,iy,p,yyy,y1,y2
        ip2 = ip2+1
        jp = ip2
        rad(jp) = dsqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
      else
        ip = ip+1
        jp = ip
      endif
      idx(np) = jp
      do  57  m = 1, 3
   57 xyz(jp,m) = xxx*xv(m) + yyy*yv(m) + z0*zv(m)
   55 continue
      if (ip2 /= np) call rx(' bug with points in hy2plt')

C --- Radial part of augmented phi's inside spheres ---
      call defrr(ophr1, np1*(lp1+1))
      call dpzero(w(ophr1),np1*(lp1+1))
      call defrr(oul,   np1*(lp1+1))
      call defrr(osl,   np1*(lp1+1))
      i1 = np-np1-np2
      if (npwr1 /= 0) then
        call us2sm2(npwr1,np1,lp1,cus1,rad(i1+1),w(oul),w(osl))
        call us2h(e1,lp1,0d0,0d0,r1,rad(i1+1),w(oul),w(osl),
     .    np1,np1,w(ophr1))
      endif
      call rlse(oul)
      call defrr(ophr2, np2*(lp2+1))
      call dpzero(w(ophr2),np2*(lp2+1))
      call defrr(oul,   np2*(lp2+1))
      call defrr(osl,   np2*(lp2+1))
      i2 = np-np2
      if (npwr2 /= 0) then
        call us2sm2(npwr2,np2,lp2,cus2,rad(i2+1),w(oul),w(osl))
        call us2h(e2,lp2,0d0,0d0,r2,rad(i2+1),w(oul),w(osl),
     .    np2,np2,w(ophr2))
      endif
      call rlse(oul)

C --- Allocate arrays for vector of phi's and xi's ---
      lmax = 0
      do  20  ix = 1, nx1
   20 lmax = max(lmax,lx1(ix))
      do  21  ix = 1, nx2
   21 lmax = max(lmax,lx2(ix))
      nlm = (lmax+1)**2
      call defrr(oxi,    np*nlm*(nx1+nx2))
      nlmp1 = (lp1+1)**2
      call defrr(oph1,   np*nlmp1)
      nlmp2 = (lp2+1)**2
      call defrr(oph2,   np*nlmp2)
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm,nlmp1,nlmp2))
C --- Tabulate xi's and phi's (Hankels) for this mesh ---
      if (iprint() >= 60) call tm('start smoothed h')
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm))
      rsmp1 = rsm1*.6d0
      rsmp2 = rsm2*.6d0
      call hy2sxi(np,xyz,nlm,
     .tau(1,1),nx1,lx1,ex1,nlmp1,rsmp1,rsm1,lp1,e1,w(ophr1),i1,np1,1,
     .tau(1,2),nx2,lx2,ex2,nlmp2,rsmp2,rsm2,lp2,e2,w(ophr2),i2,np2,1,
     .cy,w(oidx),w(owk),w(oyl),w(oph1),w(oph2),nxi,lxi,w(oxi))
      call rlse(owk)

C --- Add together fit at each point ---
      ip = 0
      do  40  ix = 1, nx
      do  40  iy = 1, ny
        ip = ip+1
        jp = idx(ip)
        do  42  m = 1, 3
        p1(m) = xyz(jp,m)-tau(m,1)
   42   p2(m) = xyz(jp,m)-tau(m,2)
        call hy2xiv(lx1,ex1,nx1,p1,rsm1,cy,lp1,e1,ph1,
     .    xi1,jp,np,nlm,0,w(oxi),nlmp1,w(oph1))
        call hy2xiv(lx2,ex2,nx2,p2,rsm2,cy,lp2,e2,ph2,
     .    xi2,jp,np,nlm,nx1,w(oxi),nlmp2,w(oph2))
        ilm1 = ilm0
        ilm2 = jlm0
        sum = 0d0
        do  31  j = 1, nf1
   31   sum = sum + c1(j,ilm1,ilm2)*xi1(j)
        do  32  j = 1, nf2
   32   sum = sum + c2(j,ilm2,ilm1)*xi2(j)
        fit(ix,iy) = sum*fac
        prd = ph1(ilm1)*ph2(ilm2)*fac
        top = dmax1(top,prd)
        bot = dmin1(bot,prd)
        prod(ix,iy) = prd
   40 continue
      write(6,561) bot,top
  561 format(' bot,top',5f12.6)

C --- Do the plot ---
C*    call begplt(25.,18.,'h')
      yszcm = 4.8
      xszcm = yszcm*(x2-x1)/(y2-y1)
      if (xszcm > 24.) then
        yszcm = yszcm*24./xszcm
        xszcm = 24.
      endif
      write(6,772) xszcm,yszcm
  772 format(' xszcm=',f12.6,'   yszcm=',f12.6)
C*    call posicm(0.1,xszcm,17.9-yszcm,17.9)
C|    call lsize(.6,.6)
      h1 = bot
      h2 = top
      dh = (h2-h1)/19.99
      h1 = -3.1d0
      h2 = +3.1d0
      dh = 0.199999
      if (dmax1(top,-bot) <= 1d0) then
        h1 = -1.05d0
        h2 = 1.05d0
        dh = 0.099999d0
      endif
      if (dmax1(top,-bot) <= 0.2d0) then
        h1 = -0.22d0
        h2 = 0.22d0
        dh = 0.03999999d0
      endif
      print *, 'h=',h1,h2,dh
      open(81,file='prod')
      open(82,file='fit')
      write(81,*) nx,ny
      write(82,*) nx,ny
      do  77  i = 1, nx
        write(81,810) (prod(i,j), j=1,ny)
        write(82,810) (fit(i,j), j=1,ny)
  810   format(5f15.7)
   77 continue
c      call pl2cp(prod(1,1),nx,ny,x1,x2,y1,y2,0d0,0d0,h1,h2,dh,cpar)
c      call pl2cp(fit(1,1),nx,ny,x1,x2,y1,y2,1.1d0,0d0,h1,h2,dh,cpar)
      stop
      end
      double precision function sp3(a,b)
      real*8 a(3),b(3)
      sp3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      end
      subroutine solhsm(r,e,a,lmax,hl,cy)
c  real smoothed solid hankel functions, negative energies
      implicit real*8 (a-h,p-z), integer(o)
      dimension r(3),hl(1),cy(1),psi(0:20)
      call sylm(r,hl,lmax,r2)
      r1 = dsqrt(r(1)**2+r(2)**2+r(3)**2)

      if (e < -1d-6) call hansmr(r1,e,a,psi,lmax)
      if (dabs(e) < 1d-6) call ropgau(1/a,lmax,1,r1,psi,000)
      if (e > 1d-6) call ropgau(e,lmax,1,r1,psi,000)
      ilm = 0
      do  10  l = 0, lmax
      nm = 2*l+1
      do  10  m = 1, nm
      ilm = ilm+1
  10  hl(ilm) = psi(l)*cy(ilm)*hl(ilm)
      end
      subroutine hy2xiv(lx,ex,nx,p,rsm,cy,lp,ep,ph,xi,
     .  ip,np,nlm,nx1,xiv,nlmp,phi)
C- Vector of strung-together xi-values from tabulated points
      implicit none
      integer lx(1),nlm,nlmp,np,ip,nx,nx1
      double precision ex(1),cy(100),hl(100),p(3),xi(1),rsm
      double precision xiv(np,nlm,1),ph(1),phi(np,nlmp),ep
      double precision diff
      integer ie,ixi,jlm,lp,lxx,i


      do  10  i = 1, nlmp
   10 ph(i) = phi(ip,i)
C     call solhsm(p,ep,1/rsm,lp,ph,cy)

      ixi = 0
      do  5  ie = 1,nx
      lxx = lx(ie)
C     call solhsm(p,ex(ie),1/rsm,lxx,hl,cy)
      do 5 jlm = 1,(lxx+1)**2
        ixi = ixi+1
        xi(ixi) = xiv(ip,jlm,nx1+ie)
C        diff = hl(jlm)-xiv(ip,jlm,nx1+ie)
C        if (dabs(diff/(hl(jlm)+1d-50)) > 1d-6)
C     .    print *, 'hy2xiv: ',ip,ie,jlm,hl(jlm), diff
C       xi(ixi) = hl(jlm)
    5 continue
      end
