      subroutine rxcrho(nxi,lxi,exi,rsm,rint,nri,rhoi,cy,n0,x,y,z,
     .    rho,rsq,q,n)
c  for one atom and xy-plane, makes rho
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),exi(1),cy(1),rhoi(nri,*),x(1),y(1),z(1),
     .  rho(n,*),q(n,2),rsq(n),ioff(10)
      real w(1)
      common /w/ w
      nsp = lsp()+1
      lmax=-1
      ioff(1)=0
      do 1 ie=1,nxi
      ioff(ie+1)=ioff(ie)+(lxi(ie)+1)**2
  1   lmax=max0(lmax,lxi(ie))
      call dpzero(rho,  n*nsp)
c ------- initialize rsq -----------------
      do 2 i=1,n
  2   rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)
c ------- make all radial parts xi -------
      call defrr(oxi,   n*(lmax+1)*nxi)
      call defrr(owk,   (4+lmax)*n)
      call defrr(oidx,  2*n)
      call hansr(rsm,0,lmax,nxi,lxi,exi,rsq,n,n,w(oidx),w(owk),
     .  000,w(oxi))
      call rlse(owk)
c ------- start big loop over m ------------
      call defrr(of,   n)
      call defrr(og,   n)
      call defrr(ocm,  n)
      call defrr(osm,  n)
      call defrr(oh,   n)
      do 10 m=0,lmax
        call ropcsm(m,n,x,y,w(oh),w(ocm),w(osm))
        do 9 isp=1,nsp
          call dpzero(w(of),  n)
          call dpzero(w(og),  n)
          do 11 l=m,lmax
            call ropqlm(m,l,n,rsq,z,q,kk)
            call ropxfg(m,l,nxi,lxi,rhoi(1,isp),cy,lmax,ioff,
     .         w(oxi),q(1,kk),w(oh),w(of),w(og),n)
  11      continue
          call roprfg(m,w(of),w(og),w(ocm),w(osm),rho(1,isp),n)
   9    continue
  10  continue
      call rlse(oxi)
      end
