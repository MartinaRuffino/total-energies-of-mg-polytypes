      subroutine pauskl(r,rsm,p,lmax,kmax,n0)
c  Make radial parts, divided by r**l, of polynomials P_kL
c  which give delta_kk'*delta_LL' when integrated times G_kL.
      implicit real*8 (a-h,p-z), integer(o)
      dimension p(0:n0,0:1)
      fpi=16d0*datan(1.d0)
      b=0.5d0*rsm*rsm
      if(kmax < 0.or.lmax < 0) return
c ---------- do explicitly for k=0,1 ------
      df=1d0
      do 7 l=0,lmax
      p(0,l)=fpi/df
      df=df*(2*l+3)
  7   if(kmax >= 1) p(1,l)=fpi*(r*r-(2*l+3)*b)/(2d0*df)
c ---------- recursion for higher k ----------
      do 8 k=2,kmax
      do 8 l=0,lmax
      z=2*k*(2*k+2*l+1)
  8   p(k,l)=(r*r*p(k-1,l)-(4*k+2*l-1)*p(k-1,l)*b-p(k-2,l)*b*b)/z
      end
