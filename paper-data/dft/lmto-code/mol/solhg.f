      subroutine solhg(e,dr,lmax,ndim,hl,ghl,cy)
C- Solid Hankel functions and gradients.
C  Note that dimension must correspond to lmax+1.
      implicit real*8 (a-h,p-z), integer (o)
      dimension cy(1),dr(1),hl(ndim),ghl(ndim,3),phi(0:30),psi(0:30)
      nlm=(lmax+1)**2
      if((lmax+2)**2 > ndim) call rx('solhg: ndim too small')

C ------ make solid Hankel functions ------
      call sylm(dr,hl,lmax+1,r2)
      call bessl(e*r2,lmax+1,phi,psi)
      ilm=0
      fac=dsqrt(r2)
      do 10 l=0,lmax+1
      fac=fac/r2
      do 10 m=-l,l
      ilm=ilm+1
  10  hl(ilm)=fac*psi(l)*cy(ilm)*hl(ilm)

C ------ make gradients ----------
      do 20 m=1,3
      do 20 ilm=1,nlm
  20  ghl(ilm,m)=0d0

      nlm1=lmax*lmax
      do 22 ilm=1,nlm
      call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
      ghl(ilm,1)=ghl(ilm,1)-cx1*hl(kx1)-cx2*hl(kx2)
      ghl(ilm,2)=ghl(ilm,2)-cy1*hl(ky1)-cy2*hl(ky2)
      ghl(ilm,3)=ghl(ilm,3)-cz*hl(kz)
      if(ilm <= nlm1) then
        ghl(kx1,1)=ghl(kx1,1)+e*cx1*hl(ilm)
        ghl(kx2,1)=ghl(kx2,1)+e*cx2*hl(ilm)
        ghl(ky1,2)=ghl(ky1,2)+e*cy1*hl(ilm)
        ghl(ky2,2)=ghl(ky2,2)+e*cy2*hl(ilm)
        ghl(kz,3)=ghl(kz,3)+e*cz*hl(ilm)
        endif
  22  continue
      end
