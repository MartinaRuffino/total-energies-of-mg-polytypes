      subroutine ftcinc(ncof,ndist,cof,ifi)
C- Read in the point-wise 2cf tables
      implicit real*8 (a-h,p-z), integer (o)
      dimension cof(ncof,ndist)
      if(iprint() >= 20) write(6,220) ncof,ndist
  220 format(/' ftcinc:  read in pointwise fits,  ncof=',i5,
     .   '   ndist=',i5)
      do 10 idist=1,ndist
  10  call dpdump(cof(1,idist),ncof,ifi)
      end
