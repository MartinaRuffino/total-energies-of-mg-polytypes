      subroutine ropqlm(m,l,n,r2,z,q,kk)
c  makes qml for m,l. must be called in sequence l=m,m+1... for fixed m
c  returns kk, which points to the current component of q.
      implicit real*8 (a-h,p-z), integer (o)
      dimension q(n,2),r2(n),z(n)
c ----- case l=m ------------------
      if(l == m) then
        a=1d0
        do 1 mm=0,m-1
  1     a=a*(2*mm+1)
        kk=1
        do 2 i=1,n
  2     q(i,kk)=a
        return
      endif
c ----- case l=m+1 ----------------
      if(l == m+1) then
        b=1d0
        do 3 mm=0,m
  3     b=b*(2*mm+1)
        kk=2
        do 4 i=1,n
  4     q(i,kk)=b*z(i)
        return
      endif
c ----- case l=m+2 and higher -----
      if(l >= m+2) then
        k2=kk
        k1=kk+1
        if(k1 == 3) k1=1
        xx=-(l+m-1d0)/(l-m)
        yy=(2*l-1d0)/(l-m)
        do 6 i=1,n
  6     q(i,k1)=xx*r2(i)*q(i,k1)+yy*z(i)*q(i,k2)
        kk=k1
        return
      endif
      end
