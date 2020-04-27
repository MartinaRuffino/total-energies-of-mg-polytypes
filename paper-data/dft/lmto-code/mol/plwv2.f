      subroutine plwv2(nbas,ib0,ips,pos,nel,rsm,nphi,lphi,ephi,n0,cy,
     .   ist,nhs,t,np,x,y,z,wfc)
c-  generate data for plot of smoothened wave function
      parameter ( nhsx=3000 )
      implicit real*8 (a-h,p-z), integer(o)
      dimension cy(1),pos(3,1),rsm(1),nphi(1),lphi(n0,1),ephi(n0,1),
     .   t(nhs,1),ips(1),wfc(np),x(np),y(np),z(np),coff(nhsx)
      real w(1)
      common /w/ w
      write(6,490) ist,np,nhs,nbas
  490 format(/'plwv2:  make wavefct for state',i5
     .       /'        np=',i5,'   nhs=',i5,'   nbas=',i4)
      if (nhs > nhsx) call rx('plwv2: increase dimen nhsx')
      call defrr(owv1,    np)
      call defrr(oq,      np*2)
      call defrr(or2,     np)

c ----- start loop over atoms -----
      do 10 ib=1,nbas
c|     if(ib0 /= 0.and.ib /= ib0) goto 10
      is=ips(ib)

c ... shift to get coordinates relative to sphere center
      write(6,470) ib,pos(1,ib),pos(2,ib),pos(3,ib),rsm(is)
  470 format('ib=',i3,'   pos=',3f8.3,'   rsm=',f8.3)
      do 11 i=1,np
      x(i)=x(i)-pos(1,ib)
      y(i)=y(i)-pos(2,ib)
  11  z(i)=z(i)-pos(3,ib)

c ... put coeffs of this atom into a vector
      ii=0
      jj=0
      do 30 je=1,nel
      do 30 jb=1,nbas
      js=ips(jb)
      nlm=(lphi(je,js)+1)**2
      do 30 jlm=1,nlm
      jj=jj+1
      if(jb == ib) then
        ii=ii+1
        coff(ii)=t(jj,ist)
        endif
  30  continue
c|    write(6,220) ib,jj,ii
c|220 format(' coffs for ib=',i5,'   at end, jj=',i4,'   ii=',i3)

c ---- make and add wfc contribution from this atom ----
      call rx('update call to rxcrho in plwv2')
C      rint=1d0
C      call rxcrho(nphi(is),lphi(1,is),ephi(1,is),rsm(is),rint,
C     .   coff,cy,n0,x,y,z,w(owv1),w(or2),w(oq),np)

      call dpadd(wfc,w(owv1),1,np,1d0)

c ... shift points back
      do 12 i=1,np
      x(i)=x(i)+pos(1,ib)
      y(i)=y(i)+pos(2,ib)
  12  z(i)=z(i)+pos(3,ib)
  10  continue

      call rlse(owv1)
      end

