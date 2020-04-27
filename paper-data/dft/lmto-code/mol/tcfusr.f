      subroutine tcfusr(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,
     .    ndim,nrhs,b,ndimx,nrhsx,cof)
c  unsort = inverse of tcfsrt; make b from cof
      implicit real*8 (a-h,p-z), integer (o)
      dimension ndim(0:20),irhs(0:20),nrhs(0:20),
     .   b(ndimx,nrhsx,0:mmax),cof(1)
      call getigv(1,ipr)
      ic=0
      do 30 m=0,mmax
  30  irhs(m)=0
      call dpzero(b,  ndimx*nrhsx*(mmax+1))
c ------- start loop over pairs phi1,phi2 ---------
      npair=0
      do 20 l1=0,lp1
      do 20 m1=0,l1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 21 l2=0,ltop
      mtop=l2
      if(lsym == 1.and.l2 == l1) mtop=m1
      do 21 m2=0,mtop
      npair=npair+1
      mp=m1+m2
      mm=iabs(m1-m2)
      fac=0.5d0
      if(m1*m2 == 0) fac=1.d0
c --------- this part for m3=m1+m2 -----------------
      if(mp <= mmax) then
        irhs(mp)=irhs(mp)+1
        jrhs=irhs(mp)
        do 31 i=1,ndim(mp)
        ic=ic+1
  31    b(i,jrhs,mp)=cof(ic)/fac
        endif
c --------- this part for m3=abs(m1-m2) ------------
      if(mm <= mmax.and.m1*m2 > 0) then
        irhs(mm)=irhs(mm)+1
        jrhs=irhs(mm)
        do 32 i=1,ndim(mm)
        ic=ic+1
  32    b(i,jrhs,mm)=cof(ic)/fac
        endif

  21  continue
  20  continue
      if(ipr >= 40) write(6,440) npair,ic
  440 format(/' tcfusr:  number of phi1*phi2 pairs=',i6
     .       /'          tot number of coefficient=',i6)

      return
      end
