      subroutine tcfxpn(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,
     .   cof,a,nxi,nlm1,nlm2)
c  expand compressed coeffs (cof) to full coeff set, a.
      implicit real*8 (a-h,p-z), integer (o)
      dimension cof(1),a(nxi,nlm1,nlm2),ipx1(20),ipx2(20),
     .   lx1(1),lx2(1)
      call getpr(ipr)
      call dpzero(a,nxi*nlm1*nlm2)
c ------ setup some pointers --------
      jpxi=0
      do 1 ie=1,nx1
      ipx1(ie)=jpxi
  1   jpxi=jpxi+(lx1(ie)+1)**2
      do 2 ie=1,nx2
      ipx2(ie)=jpxi
  2   jpxi=jpxi+(lx2(ie)+1)**2
      ic=0
c ------- start loop over pairs phi1,phi2 ---------
      do 20 l1=0,lp1
      do 20 m1=0,l1
      ilm10=l1*(l1+1)+1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 21 l2=0,ltop
      mtop=l2
      if(lsym == 1.and.l2 == l1) mtop=m1
      do 21 m2=0,mtop
      ilm20=l2*(l2+1)+1
      sig=1.d0
      if(m2 > m1) sig=-1.d0
c --------- this part for m3=m1+m2, site 1 ----------
      m3=m1+m2
      if(m3 <= mmax) then
        do 31 ie=1,nx1
        do 31 l3=m3,lx1(ie)
        ioff=ipx1(ie)+l3*(l3+1)+1
        ic=ic+1
                       a(ioff+m3,ilm10+m1,ilm20+m2)= cof(ic)
        if(m1 > 0)    a(ioff-m3,ilm10-m1,ilm20+m2)= cof(ic)
        if(m2 > 0)    a(ioff-m3,ilm10+m1,ilm20-m2)= cof(ic)
        if(m1*m2 > 0) a(ioff+m3,ilm10-m1,ilm20-m2)=-cof(ic)
        if(lsym == 1) then
        isl=1-2*mod(l1+l2+l3,2)
        joff=ipx2(ie)+l3*(l3+1)+1
                       a(joff+m3,ilm20+m2,ilm10+m1)= cof(ic)*isl
        if(m1 > 0)    a(joff-m3,ilm20+m2,ilm10-m1)= cof(ic)*isl
        if(m2 > 0)    a(joff-m3,ilm20-m2,ilm10+m1)= cof(ic)*isl
        if(m1*m2 > 0) a(joff+m3,ilm20-m2,ilm10-m1)=-cof(ic)*isl
        endif
  31    continue
c --------- this part for m3=m1+m2, site 2 ----------
        do 32 ie=1,nx2
        do 32 l3=m3,lx2(ie)
        ioff=ipx2(ie)+l3*(l3+1)+1
        ic=ic+1
                       a(ioff+m3,ilm10+m1,ilm20+m2)= cof(ic)
        if(m1 > 0)    a(ioff-m3,ilm10-m1,ilm20+m2)= cof(ic)
        if(m2 > 0)    a(ioff-m3,ilm10+m1,ilm20-m2)= cof(ic)
        if(m1*m2 > 0) a(ioff+m3,ilm10-m1,ilm20-m2)=-cof(ic)
        if(lsym == 1) then
        isl=1-2*mod(l1+l2+l3,2)
        joff=ipx1(ie)+l3*(l3+1)+1
                       a(joff+m3,ilm20+m2,ilm10+m1)= cof(ic)*isl
        if(m1 > 0)    a(joff-m3,ilm20+m2,ilm10-m1)= cof(ic)*isl
        if(m2 > 0)    a(joff-m3,ilm20-m2,ilm10+m1)= cof(ic)*isl
        if(m1*m2 > 0) a(joff+m3,ilm20-m2,ilm10-m1)=-cof(ic)*isl
        endif
  32    continue
      endif
c --------- this part for m3=abs(m1-m2), site 1 ------------
      m3=iabs(m1-m2)
      if(m3 <= mmax.and.m1*m2 > 0) then
        do 33 ie=1,nx1
        do 33 l3=m3,lx1(ie)
        ioff=ipx1(ie)+l3*(l3+1)+1
        ic=ic+1
                    a(ioff+m3,ilm10+m1,ilm20+m2)= cof(ic)
        if(m3 > 0) a(ioff-m3,ilm10-m1,ilm20+m2)= cof(ic)*sig
        if(m3 > 0) a(ioff-m3,ilm10+m1,ilm20-m2)=-cof(ic)*sig
                    a(ioff+m3,ilm10-m1,ilm20-m2)= cof(ic)
        if(lsym == 1) then
        isl=1-2*mod(l1+l2+l3,2)
        joff=ipx2(ie)+l3*(l3+1)+1
                    a(joff+m3,ilm20+m2,ilm10+m1)= cof(ic)*isl
        if(m3 > 0) a(joff-m3,ilm20+m2,ilm10-m1)= cof(ic)*sig*isl
        if(m3 > 0) a(joff-m3,ilm20-m2,ilm10+m1)=-cof(ic)*sig*isl
                    a(joff+m3,ilm20-m2,ilm10-m1)= cof(ic)*isl
        endif
  33    continue
c --------- this part for m3=abs(m1-m2), site 2 ------------
        do 34 ie=1,nx2
        do 34 l3=m3,lx2(ie)
        ioff=ipx2(ie)+l3*(l3+1)+1
        ic=ic+1
                    a(ioff+m3,ilm10+m1,ilm20+m2)= cof(ic)
        if(m3 > 0) a(ioff-m3,ilm10-m1,ilm20+m2)= cof(ic)*sig
        if(m3 > 0) a(ioff-m3,ilm10+m1,ilm20-m2)=-cof(ic)*sig
                    a(ioff+m3,ilm10-m1,ilm20-m2)= cof(ic)
        if(lsym == 1) then
        isl=1-2*mod(l1+l2+l3,2)
        joff=ipx1(ie)+l3*(l3+1)+1
                    a(joff+m3,ilm20+m2,ilm10+m1)= cof(ic)*isl
        if(m3 > 0) a(joff-m3,ilm20+m2,ilm10-m1)= cof(ic)*sig*isl
        if(m3 > 0) a(joff-m3,ilm20-m2,ilm10+m1)=-cof(ic)*sig*isl
                    a(joff+m3,ilm20-m2,ilm10-m1)= cof(ic)*isl
        endif
  34    continue
      endif

  21  continue
  20  continue
c ---------------------------------------
c|99  ilm1=0
c|    write(6,*) 'enter ilm1,ilm2'
c|    read (5,*)  ilm1,ilm2
c|    if(ilm1 <= 0) return
c|    ixi=0
c|    do 5 ie=1,nx1
c|    do 5 jlm=1,(lx1(ie)+1)**2
c|    ixi=ixi+1
c|    aa=a(ixi,ilm1,ilm2)
c|5   if(dabs(aa) > 1.d-6) write(6,800) ixi,1,ie,jlm,aa
c|    do 6 ie=1,nx2
c|    do 6 jlm=1,(lx2(ie)+1)**2
c|    ixi=ixi+1
c|    aa=a(ixi,ilm1,ilm2)
c|6   if(dabs(aa) > 1.d-6) write(6,800) ixi,2,ie,jlm,aa
c|800 format(i5,3x,3i6,f14.5)
c|    if(.true.) goto 99

      return
      end
