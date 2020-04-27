      subroutine parsgn(t,nt,g,a,ng)
c  parses the string t containing the group element describers
c  and makes a (3x3) matrix g and vector a for each element.
c  space group ops are defined as:  (g,a)(p)=g*p+a
c
c  point group element describers:
c    rnv = n-fold rot around vec v, n pos integer   ) also products
c    mv  = mirror taking vec v into -v              ) of form 'a*b' of
c    i   = inversion                                ) any two of these
c
c  glide planes, screw axes etc:
c    append  ':v'  to any point operation above where v is the
c    translation vector for the operation.
c
c  vectors v are all given as (x,y,z) with real numbers x,y,z
c  or abbreviate   'd' for (1,1,1), 'x' for (1,0,0),  y,z  similar.
c
      character*1 t(nt)
      real*8 g(9,1),h(9),a(3,1)
      ng=0
      i=1
  90  call parsbl(t,nt,i)
      if(i > nt) return
      ng=ng+1
      call parsop(t,nt,i,g(1,ng))
      if(t(i) == '*') then
        i=i+1
        call parsop(t,nt,i,h)
        call grpprd(g(1,ng),h,g(1,ng))
        endif
      do 1 m=1,3
  1   a(m,ng)=0.d0
      if(t(i) == ':') then
        i=i+1
        call parsvc(t,nt,i,a(1,ng))
        endif
      if(t(i) /= ' ') call rx('parsgn: expect blank or ":" after op')
      goto 90
      end
c ------------------------
      subroutine parsop(t,nt,i,a)
      real*8 v(3),sp,c,s,pi,a(3,3)
      character*1 t(nt)
      pi=4.d0*datan(1.d0)
      if(t(i) == 'r') then
        i=i+1
        call parspi(t,nt,i,nrot,ndig)
        if(ndig == 0) call rx('parsop: expect integer after "r"')
        call parsvc(t,nt,i,v)
        sp=0.d0
        do 13 k=1,3
  13    sp=sp+v(k)*v(k)
        sp=1.d0/dsqrt(sp)
        do 14 k=1,3
  14    v(k)=v(k)*sp
        c=dcos(2.d0*pi/nrot)
        s=dsin(2.d0*pi/nrot)
        do 16 k=1,3
        do 15 j=1,3
  15    a(j,k)=(1.d0-c)*v(j)*v(k)
  16    a(k,k)=a(k,k)+c
        a(1,2)=a(1,2)+s*v(3)
        a(3,1)=a(3,1)+s*v(2)
        a(2,3)=a(2,3)+s*v(1)
        a(2,1)=a(2,1)-s*v(3)
        a(1,3)=a(1,3)-s*v(2)
        a(3,2)=a(3,2)-s*v(1)
      else if(t(i) == 'm') then
        i=i+1
        call parsvc(t,nt,i,v)
        sp=0.d0
        do 10 k=1,3
  10    sp=sp+v(k)*v(k)
        do 11 j=1,3
        do 12 k=1,3
  12    a(k,j)=-2.d0*v(k)*v(j)/sp
  11    a(j,j)=a(j,j)+1.d0
      else if(t(i) == 'i') then
        do 35 k=1,3
        do 34 j=1,3
  34    a(k,j)=0.d0
  35    a(k,k)=-1.d0
        i=i+1
      else
        call rx('parsop: op must start with "r","m" or "i"')
      endif
      return
      end
c ------------------------------------
      subroutine parsvc(t,nt,i,v)
      real*8 r,v(3)
      character*1 t(nt)
      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0
      if(t(i) == 'x') v(1)=1.d0
      if(t(i) == 'y') v(2)=1.d0
      if(t(i) == 'z') v(3)=1.d0
      if(t(i) == 'x'.or.t(i) == 'y'.or.t(i) == 'z') goto 90
      if(t(i) == 'd') then
        v(1)=1.d0
        v(2)=1.d0
        v(3)=1.d0
        goto 90
        endif
      if(t(i) /= '(') call rx('parsvc: wrong vector format')
      do 10 m=1,3
      i=i+1
      call parssr(t,nt,i,r)
      if(m < 3.and.t(i) /= ',') call rx('parsvc: "," expected')
  10  v(m)=r
      if(t(i) /= ')') call rx('parsvc: ")" expected')
  90  i=i+1
      return
      end
c ----------------------------------
      subroutine parspi(t,nt,i,integr,ndig)
      character*1 t(nt),csym(10)
      data csym/'0','1','2','3','4','5','6','7','8','9'/
      integr=0
      ndig=0
      i=i-1
  99  i=i+1
      do 1 is=1,10
      if(csym(is) == t(i)) then
        integr=10*integr+is-1
        ndig=ndig+1
        goto 99
        endif
  1   continue
      return
      end
c ----------------------------------
      subroutine parssr(t,nt,i,r)
      character*1 t(nt),csym(10)
      real*8 r
      isig=1
      if(t(i) == '-') isig=-1
      if(t(i) == '-'.or.t(i) == '+') i=i+1
      call parspi(t,nt,i,intg1,ndig)
      if(t(i) == '.') then
        i=i+1
        call parspi(t,nt,i,intg2,ndig)
        r=isig*(intg1+intg2*(0.1d0**ndig))
      else
        r=intg1*isig
      endif
      return
      end
c ----------------------------------
      subroutine parsbl(t,nt,i)
      character*1 t(nt)
  99  if(t(i) /= ' ') return
      i=i+1
      if(i > nt) return
      goto 99
      end
c --------------------------
      subroutine grpprd(g1,g2,g1xg2)
      real*8 g1(3,3),g2(3,3),g1xg2(3,3),h(3,3),sum
      do 10 i=1,3
      do 10 j=1,3
      sum=0.d0
      do 11 k=1,3
  11  sum=sum+g1(i,k)*g2(k,j)
  10  h(i,j)=sum
      do 12 j=1,3
      do 12 i=1,3
  12  g1xg2(i,j)=h(i,j)
      return
      end
