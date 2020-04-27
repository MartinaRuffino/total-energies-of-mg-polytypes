      subroutine tcfrt1(d,rmat,lmax,ndim)
c  makes matrix to rotate tcf to a general direction
      implicit real*8 (a-h,o-z)
      dimension d(3),rmat(ndim,ndim)
      nlm=(lmax+1)**2
      if(nlm > ndim) call rx('tcfrt1: nlm > ndim''')
      call getpr(ipr)
c  ---- alfa,beta,gama=euler angles to rotate d into z-axis -----
      pi=4d0*datan(1d0)
      d1=dsqrt(d(1)**2+d(2)**2+d(3)**2)
      d2=dsqrt(d(1)**2+d(2)**2)
C|       write(6,222) d,d1,d2,d2/d1
C|222    format(' d',3f10.5,'  d1,d2'2f10.5,'  d2/d1',1p,d12.3)
      if(d2/d1 > 1d-7) then
        beta=-dacos(d(3)/d1)
        alfa=-datan2(d(2),d(1))
      else
        alfa=0.d0
        if(d(3) > 0d0) beta=0d0
        if(d(3) < 0d0) beta=-pi
      endif
      gama=0.d0
      if(ipr > 80) write(6,873) d,alfa,beta,lmax
C|                  write(6,873) d,alfa,beta,lmax
  873 format(' tcfrt1:  d=',3f8.3,'   alf,bet=',2f8.3,'   lmax=',i2)
      do 7 jlm=1,nlm
      do 7 ilm=1,nlm
  7   rmat(ilm,jlm)=0.d0
      do 8 l=0,lmax
      ilm1=l*l+1
c?8   call ylmrot(l,alfa,beta,gama,rmat(ilm1,ilm1),nlm)
  8   call ylmrot(l,alfa,beta,gama,rmat(ilm1,ilm1),ndim)

      return
      end
