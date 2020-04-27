      subroutine rxcrh1(nxi,lxi,exi,rsm,rint,nri,rhoi,cg,jcg,indxcg,
     .    nrx,x,y,z,rho,grho,g2rho,ggrho,n)
C- for one atom and xy-plane, makes rho
Cu   21 Mar 07 (S. Lozovi) First created
      implicit real*8 (a-h,p-z), integer (o)
      dimension lxi(1),exi(1),rhoi(nri,*),x(1),y(1),z(1),
     .  cg(1),jcg(1),indxcg(1),
     .  rho(nrx,*),grho(nrx,3,*),g2rho(nrx,5,*),ggrho(nrx,*),
     .  ioff(10)
c  Local variables
      dimension lxig(10)
      real w(1)
      common /w/ w

      call tcn('rxcrh1')

c --------- Sanity check ---------
      if(nxi >= 10) call rx('rxcrh1: increase dimension of ioff')

c --------- If GGA we need kmax=2 and lmax -> lmax+2  ---------
      if(lxcg() == 0) then
       kmax=0
      else
       kmax=2
      endif

c --------- Find lmax and make the gga-copy of lxi ---------
      lmax = -10
      ioff(1) = 0
      do  ie = 1,nxi
        lxig(ie)=lxi(ie)+kmax
        ioff(ie+1) = ioff(ie)+(lxi(ie)+1)**2
        lmax = max0(lmax,lxi(ie))
      enddo
      lmaxg=lmax+kmax

c ------- harmonic polynomials for the whole plane -------
      call defrr(oyl,     (lmaxg+1)**2*nrx)
      call defrr(or2,     n)
      call ropyln(n,x,y,z,lmaxg,nrx,w(oyl),w(or2))
c this is now done inside ropyln
c ------- initialize rsq -----------------
c     do 2 i=1,n
c 2   rsq(i)=x(i)*x(i)+y(i)*y(i)+z(i)*z(i)

c ---- if GGA, make energy-independent Gaussians -----
      if(kmax > 0) then
        call defrr(og,     (lmaxg+1)*(kmax+1)*nrx)
        call gklv(rsm,w(or2),nrx,kmax-1,n,lmaxg,kmax-1,w(og))
      else
        og = 1
      endif

c ------- make radial parts of smooth Hankels, xi -------
      call defrr(oxi,   nrx*(lmaxg+1)*nxi)
      call defrr(owk,   (4+lmaxg)*nrx)
      call defrr(oidx,  2*nrx)
      call hansr(rsm,0,lmaxg,nxi,lxig,exi,w(or2),nrx,n,w(oidx),w(owk),
     .  000,w(oxi))
      call rlse(owk)

c -- make solid hankels and, if GGA, their gradients, Laplacians, and second derivatives --
        call defrr(oh,     (lmax+1)**2*nrx*nxi)
      if(kmax > 0) then
        call defrr(ogh,   -(lmax+1)**2*nrx*3*nxi)
        call defrr(oggh,   (lmax+1)**2*nrx*nxi)
        call defrr(og2h,  -(lmax+1)**2*nrx*5*nxi)
      else
        ogh = 1
        oggh = 1
        og2h = 1
      endif
c allocate a work array to keep hankels up to lmax+2
      call defrr(ohh,    (kmax+1)*(lmaxg+1)**2*nrx)

      call hdercg(rsm,nrx,nxi,lmaxg,lxig,lmax,lxi,exi,w(oxi),w(oyl),
     .  kmax,cg,jcg,indxcg,n,w(og),w(ohh),w(oh),w(ogh),w(og2h),w(oggh))

      call rlse(ohh)

c ------ make electronic density ---------
      call rhosum(nrx,nxi,lxi,nri,rhoi,lmax,ioff,
     .   w(oh),w(ogh),w(og2h),w(oggh),n,rho,grho,g2rho,ggrho)

c ------- old code
cc ------- start big loop over m ------------
c      call defrr(of,   n)
c      call defrr(og,   n)
c      call defrr(ocm,  n)
c      call defrr(osm,  n)
c      call defrr(oh,   n)
c      do 10 m=0,lmax
c        call ropcsm(m,n,x,y,w(oh),w(ocm),w(osm))
c        do 9 isp=1,nsp
c          call dpzero(w(of),  n)
c          call dpzero(w(og),  n)
c          do 11 l=m,lmax
cc           call ropqlm(m,l,n,rsq,z,q,kk)
c            call ropqlm(m,l,n,w(or2),z,q,kk)
c            call ropxfg(m,l,nxi,lxi,rhoi(1,isp),cy,lmax,ioff,
c     .         w(oxi),q(1,kk),w(oh),w(of),w(og),n)
c  11      continue
c          call roprfg(m,w(of),w(og),w(ocm),w(osm),rho(1,isp),n)
c   9    continue
c  10  continue
cc     call rlse(oxi)

      call rlse(oyl)
      call tcx('rxcrh1')
      end
