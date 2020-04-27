      subroutine slhhpg(e,lmax,ndim,hl,ghl,hd,ghd)
C- Solid Hankel functions, energy derivatives and gradients given hl,hd
      implicit real*8 (a-h,p-z), integer (o)
      dimension hl(1),ghl(ndim,3),hd(1),ghd(ndim,3)
      nlm=(lmax+1)**2
      if((lmax+2)**2 > ndim) call rx('solhgp: ndim too small')

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
