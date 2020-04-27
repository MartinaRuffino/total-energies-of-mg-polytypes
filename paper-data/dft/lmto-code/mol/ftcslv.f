      subroutine ftcslv(mmax,s,b,ndim,nrhs,ndimx,nsdmx,nrhsx)
C- Solve lsq-fit problem for tcf. output: coeffs in b.
      implicit real*8 (a-h,p-z), integer (o)
      dimension s(nsdmx,0:mmax),b(ndimx,nrhsx,0:mmax),
     .   ndim(0:20),nrhs(0:mmax)
      real w(1)
      common /w/ w
      if(iprint() >= 30) write(6,500) mmax,(ndim(m),m=0,mmax)
  500 format(/' ftcslv:  mmax=',i3,'    ndim=',8i5)
      call defrr(ow,   ndimx)

c --------- start loop over m -------------------
      do 10 m=0,mmax
      nd=ndim(m)
      call dspfa(s(1,m),nd,w(ow),ierr)
      if (ierr /= 0)
     .  call rx('ftcslv: matrix not pos def for m='//char(m+ichar('0')))
      do  11  jrhs = 1, nrhs(m)
        call dspsl(s(1,m),nd,w(ow),b(1,jrhs,m))
      if(iprint() >= 80) write(6,220) jrhs,(b(i,jrhs,m),i=1,nd)
  220 format(' jrhs=',i5/(1x,8f10.5))
  11  continue
  10  continue
      call rlse(ow)
      end
