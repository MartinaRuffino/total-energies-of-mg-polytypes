      subroutine ftccmp(ncof,ndist,cuf,cof,erravg,errmax,i0)
c  compare two cof arrays
      implicit real*8 (a-h,p-z), integer (o)
      dimension cuf(ncof,ndist),cof(ncof,ndist)
      call getpr(ipr)
      errmax=0d0
      erravg=0d0
      i0=0
      do 10 i=1,ncof
      top=0d0
      dif=0d0
      do 11 id=1,ndist
      top=dmax1(top,dabs(cof(i,id)))
  11  dif=dmax1(dif,dabs(cuf(i,id)-cof(i,id)))
      relerr=100d0*dif/top
      erravg=erravg+relerr/ncof
      if(relerr >= errmax) then
        if(ipr >= 35) write(6,837) i,top,dif,relerr
  837   format(10x,i5,'   top,dif',2f12.6,'   (',f6.3,' % )')
        errmax=relerr
        i0=i
        endif
  10  continue
      end
