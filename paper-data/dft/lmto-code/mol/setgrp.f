      subroutine setgrp(symgrp,ng,og,oag)
      implicit real*8 (a-h,p-z), integer (o)
      dimension qb(3,3)
      character*1 symgrp(60)
      real w(1)
      common /w/ w
      data qb /0.01d0,0d0,0d0, 0d0,0.01d0,0d0, 0d0,0d0,0.01d0/
      if(iprint() >= 20) write(6,100) symgrp
  100 format(/' setgrp: ',60a1)
      ngmx=50
      ngnmx=10
      call defrr(og,     9*ngmx)
      call defrr(oag,    3*ngmx)
      call defrr(ogn,   9*ngnmx)
      call defrr(oagn,  3*ngnmx)
      call parsgn(symgrp,60,w(ogn),w(oagn),ngen)
      call sgroup(1,w(ogn),w(oagn),ngen,w(og),w(oag),ng,ngmx,qb)
      if(iprint() >= 20) write(6,650) ng
  650 format(' number of elements in space group=',i4)
      call rlse(ogn)
CL    write(71,710) ng,symgrp
  710 format(' mce  ng',i4,1x,60a1)
      end
