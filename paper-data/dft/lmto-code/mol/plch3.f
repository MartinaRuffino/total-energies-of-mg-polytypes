      subroutine plch3(rmt,a,nr,lmxl,nbas,pos,ips,cy,orho,
     .   np,x,y,z,rho)
c- puts the numeric sphere density into array rho for points in mt-sph.
      implicit real*8 (a-h,p-z), integer (o)
      dimension ips(1),pos(3,1),cy(1),rmt(1),lmxl(1),a(1),nr(1),
     .  x(np),y(np),z(np),rho(np),p0(3),orho(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('plch3')
      write(6,983) np
  983 format(' plch3:  np=',i5)
c ------ start loop over atoms -------------------
      do 10 ib=1,nbas
      is=ips(ib)
      oro=orho(ib)
      lmaxl=lmxl(is)
      np1=0
      do 11 ip=1,np
      p0(1)=x(ip)-pos(1,ib)
      p0(2)=y(ip)-pos(2,ib)
      p0(3)=z(ip)-pos(3,ib)
      d1=dsqrt(p0(1)**2+p0(2)**2+p0(3)**2)
      if(d1 < rmt(is)) then
        call yypch3(w(oro),lmaxl,a(is),rmt(is),nr(is),cy,p0,rhop)
        rho(ip)=rhop
        np1=np1+1
        endif
  11  continue
      if(np1 > 0.and.ipr >= 20) write(6,720) ib,np1
  720 format(' ib=',i5,'   np1=',i5)
  10  continue
      call tcx('plch3')
      end

c ------ sub yypch3 --------------
      subroutine yypch3(rhol,lmax,a,rmt,nr,cy,p0,rho)
c  evaluate 1-center expansion of rho inside sphere
      implicit real*8 (a-h,p-z), integer(o)
      dimension rhol(nr,1),cy(1),p(3),yl(300),p0(3)
      nlm=(lmax+1)**2
      b=rmt/(dexp(a*nr-a)-1.d0)
      pp=dsqrt(p0(1)**2+p0(2)**2+p0(3)**2)+1.d-6
      do 4 k=1,3
  4   p(k)=p0(k)/pp
      call sylm(p,yl,lmax,pp1)
      sum1=0.d0
      do 12 ilm=1,nlm
      call pppf2(pp,a,b,nr,rhol(1,ilm),val)
  12  sum1=sum1+cy(ilm)*yl(ilm)*val
      rho=sum1/(pp*pp)
      return
      end
c ------ pppf2 --------------------
      subroutine pppf2(r,a,b,n,f,val)
c  interpolates f on rad mesh
      implicit real*8 (a-h,p-z), integer(o)
      dimension f(1),gg(0:4)
      ff=dlog(r/b+1.d0)/a+1.d0
      n0=ff-2.d0
      n0=max0(n0,1)
      n0=min0(n0,n-4)
      x=ff-n0
      do 1 j=0,4
  1   gg(j)=f(n0+j)
      call polft4(gg,x,vul,slu,sum)
      val=vul
      return
      end
