      subroutine hyfipt(d,d0,adec,nalf,nph,mmax,
     .  ndimx,nrhsx,ndim,nrhsj,nx1,ba,b)
C- Chebyshev interpolation of HYF coefficients for sets of phi's
      implicit none
      integer mmax,ndim(0:mmax,0:1),nrhsj(0:mmax,1),nalf,nx1,
     .  ndimx,nrhsx,nph
      double precision adec,x,y,
     .  b(ndimx,nrhsx,0:mmax),ba(ndimx,nrhsx,0:mmax,nalf),
     .  d,d0,f(40)
      integer ialf,id,ir,m,nr,iph,nd1,nd2,iprint

      x = dexp(adec*(d0-d))
      y = 2*x-1
      if(iprint() >= 30) write(6,440) d,d0,adec,x,nalf
  440 format(/' hyfipt:  d=',f9.6,'   d0=',f9.6,'   adec=',f7.4,
     .   '   x=',f9.5,'   nalf=',i2)
      call dpzero(b, ndimx*nrhsx*(mmax+1))

C --- Interpolate coffs for each set of phis ---
      do  10  iph = 1, nph

C ...  Set up interpolating functions f
        do  20  ialf = 1, nalf
   20   f(ialf) = dcos((ialf-1)*dacos(y))

C ...   Combine alfa-coffs for coffs at this distance
        do  40  m = 0, mmax
        nd1 = ndim(m,nx1)
        nd2 = ndim(m,0)
        nr = nrhsj(m,iph+1)
        do  40  ir = nrhsj(m,iph)+1, nr
        do  40  ialf = 1, nalf
          do  42  id = 1, nd1
   42     b(id,ir,m) = b(id,ir,m) + f(ialf)*ba(id,ir,m,ialf)
C          if (ir == 1 .and. m == 0)
C     .      print *, 'HYFIPT: f,ba(1)  ', f(ialf),ba(1,ir,m,ialf)
          do  43  id = nd1+1, nd2
   43     b(id,ir,m) = b(id,ir,m) + f(ialf)*ba(id,ir,m,ialf)
C          if (ir == 1 .and. m == 0)
C     .      print *, 'HYFIPT: f,ba(nd2)', f(ialf),ba(nd2,ir,m,ialf)
   40   continue
   10 continue

C      print *,'HYFIPT: b(1),b(nd2)',ndim(0,0),b(1,1,0),b(ndim(0,0),1,0)

      end
