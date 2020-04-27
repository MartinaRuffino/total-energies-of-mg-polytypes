      subroutine rophs0(e,rsm,lmax,n,r,xi,job)
c  makes vectors of smoothed hankel functions.
c  this sub makes hankels by calling hansmr - non-vectorizable.
c  set job=1 to multiply xi_l by r**l.
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(n),xi(n,0:lmax),f(0:20),phi(0:20)
c ------- here if smoothed hankels ---------
      if(rsm > 1.d-1) then
        a=1d0/rsm
        do 10 i=1,n
        call hansmr(r(i),e,a,f,lmax)
        do 11 l=0,lmax
  11    xi(i,l)=f(l)
  10    continue
c ------- here if normal hankels ---------
      else
        do 20 i=1,n
        rr=1d0
        if(r(i) > 1d-20) rr=r(i)
        call bessl(e*rr**2,lmax,phi,f)
        do 21 l=0,lmax
  21    xi(i,l)=f(l)/rr**(2*l+1)
  20    continue
      endif
c ------- scale if job=1 -------------
      if(job == 1) then
        do 3 l=1,lmax
        do 3 i=1,n
  3     xi(i,l)=xi(i,l)*r(i)**l
      endif
      end
