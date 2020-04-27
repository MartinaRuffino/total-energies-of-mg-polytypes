      subroutine solhpg(e,dr,lmax,ndim,hl,ghl,hd,ghd,cy)
C- Solid Hankel functions with energy derivatives and gradients
      implicit real*8 (a-h,p-z), integer (o)
      dimension cy(*),dr(3),hl(ndim),ghl(ndim,3),phi(0:30),psi(0:30),
     .   hd(ndim),ghd(ndim,3)
      nlm=(lmax+1)**2
      if((lmax+2)**2 > ndim) call rx('solhpg: ndim too small')

C --- Make solid Hankel functions HL ---
      call sylm(dr,hl,lmax+1,r2)
      call bessl(e*r2,lmax+2,phi,psi)
      ilm=0
      fac=dsqrt(r2)
      do 10 l=0,lmax+1
        fac=fac/r2
        psidot=((l+l+1)*psi(l)-psi(l+1))/(e+e)
        do 10 m=-l, l
        ilm=ilm+1
        hd(ilm)=fac*psidot*cy(ilm)*hl(ilm)
  10    hl(ilm)=fac*psi(l)*cy(ilm)*hl(ilm)

C ------ make gradients ----------
      do 20 m=1,3
      do 20 ilm=1,nlm
      ghd(ilm,m)=0d0
  20  ghl(ilm,m)=0d0

      nlm1=lmax*lmax
      do 22 ilm=1,nlm
      call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
      ghl(ilm,1)=ghl(ilm,1)-cx1*hl(kx1)-cx2*hl(kx2)
      ghl(ilm,2)=ghl(ilm,2)-cy1*hl(ky1)-cy2*hl(ky2)
      ghl(ilm,3)=ghl(ilm,3)-cz*hl(kz)
      ghd(ilm,1)=ghd(ilm,1)-cx1*hd(kx1)-cx2*hd(kx2)
      ghd(ilm,2)=ghd(ilm,2)-cy1*hd(ky1)-cy2*hd(ky2)
      ghd(ilm,3)=ghd(ilm,3)-cz*hd(kz)
      if(ilm <= nlm1) then
        xx=e*hl(ilm)
        ghl(kx1,1)=ghl(kx1,1)+cx1*xx
        ghl(kx2,1)=ghl(kx2,1)+cx2*xx
        ghl(ky1,2)=ghl(ky1,2)+cy1*xx
        ghl(ky2,2)=ghl(ky2,2)+cy2*xx
        ghl(kz,3)=ghl(kz,3)+cz*xx
        xx=hl(ilm)+e*hd(ilm)
        ghd(kx1,1)=ghd(kx1,1)+cx1*xx
        ghd(kx2,1)=ghd(kx2,1)+cx2*xx
        ghd(ky1,2)=ghd(ky1,2)+cy1*xx
        ghd(ky2,2)=ghd(ky2,2)+cy2*xx
        ghd(kz,3)=ghd(kz,3)+cz*xx
        endif
  22  continue
      end
