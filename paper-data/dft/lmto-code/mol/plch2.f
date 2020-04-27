      subroutine plch2(nbas,ips,pos,rsm,nxi,lxi,exi,n0,ioff,cy,
     .   rhoi,np,x,y,z,rho)
c-  adds smooth interstitial rho for points in list
      implicit real*8 (a-h,p-z), integer(o)
      dimension cy(1),pos(3,1),rsm(1),nxi(1),lxi(n0,1),exi(n0,1),
     .   rhoi(1),ips(1),rho(np),x(np),y(np),z(np),ioff(1)
      real w(1)
      common /w/ w
      write(6,928) np
  928 format(' plch2:  np=',i5)
      call defrr(orho1,   np)
      call defrr(oq,      np*2)
      call defrr(or2,     np)
      do 10 ib=1,nbas
      is=ips(ib)
c ... shift to get coordinates relative to sphere center
      do 11 i=1,np
      x(i)=x(i)-pos(1,ib)
      y(i)=y(i)-pos(2,ib)
  11  z(i)=z(i)-pos(3,ib)
c ... make and add rho of this atom
      ir=ioff(ib)+1
      call rx('update call to rxcrho in plwv2')

C      rint=1d0
Cc|       call wkchk('before rxcrho')
C      call rxcrho(nxi(is),lxi(1,is),exi(1,is),rsm(is),rint,
C     .   rhoi(ir),cy,n0,x,y,z,w(orho1),w(or2),w(oq),np)

      call dpadd(rho,w(orho1),1,np,1d0)
c ... shift points back
      do 12 i=1,np
      x(i)=x(i)+pos(1,ib)
      y(i)=y(i)+pos(2,ib)
  12  z(i)=z(i)+pos(3,ib)
  10  continue
      call rlse(orho1)
      end
