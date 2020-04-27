      subroutine gauxpn(taub,taua,rsmb,rsmba,nlmb,nlma,npwr,
     .   cg,indxcg,jcg,cy,c)
c  Makes the coefficients to expand the gaussian to rsmb at taub
c  into polynomials P_kL to rsma at point taua.
      parameter( lx=16, lpx=52 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension c(0:npwr,nlma,nlmb),cy(1),cg(1),indxcg(1),jcg(1),dr(3),
     .  yl((lx+1)**2),phi(0:lx),taua(3),taub(3),g(0:lpx,0:lx)
      fpi=16.d0*datan(1.d0)
      gamh=0.25d0*rsmb*rsmb
      gamg=0.25d0*rsmba*rsmba
      gam=gamh+gamg
      a=0.5d0/dsqrt(gam)
      do 1 m=1,3
    1 dr(m)=taub(m)-taua(m)
      lmaxa=ll(nlma)
      lmaxb=ll(nlmb)
      lmmax=lmaxa+lmaxb
      kmax=max0(lmaxa,lmaxb)+npwr
      if(lmmax > lx) call rx('gauxpn: increase lx')
      if(kmax > lpx) call rx('gauxpn: increase lpx')
      call sylm(dr,yl,lmmax,r2)
      do 2 ilm=1,(lmmax+1)**2
    2 yl(ilm)=yl(ilm)*cy(ilm)
      r=dsqrt(r2)
C --- Set up g(k,l) --------------
      call gauskl(r,0d0,a,g,lmmax,kmax,lpx)
C --- Combine with Clebsch-Gordan coefficients ---
      call dpzero(c,(npwr+1)*nlmb*nlma)
      do 11 ila=1,nlma
      la=ll(ila)
      do 11 ilb=1,nlmb
      lb=ll(ilb)
      ii=max0(ila,ilb)
      indx=(ii*(ii-1))/2+min0(ila,ilb)
      icg1=indxcg(indx)
      icg2=indxcg(indx+1)-1
      do 11 icg=icg1,icg2
      ilm=jcg(icg)
      lm=ll(ilm)
      k=(la+lb-lm)/2
      fac=(-1d0)**lb*cg(icg)*yl(ilm)
      do 11 ip=0,npwr
   11 c(ip,ila,ilb)=c(ip,ila,ilb)+fac*g(k+ip,lm)

      end
