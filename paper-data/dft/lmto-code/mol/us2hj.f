      subroutine us2hj(e,lmax,el,tl,rmt,rofi,ul,sl,nrx,nr,job,xi,b)
C- Linear combinations of (u,s) that match differentiably onto e
C  xi and b: true radial parts/r**l so the solid hankel functions are
C  hl(ilm) = xi(l)*cy(ilm)*yl(ilm) or jl(ilm) = b(l)*cy(ilm)*yl(ilm)
C  job=1: true radial parts for xi and b
      implicit none
      integer nrx,nr,lmax,job
      double precision xi(nrx,0:1),b(nrx,0:1),ul(nr,0:1),sl(nr,0:1),
     .  rofi(1),el,tl,e,rmt
      integer ir,l,i1,il
      double precision phi(0:20),psi(0:20),sh,sj,vh,vj,rl,er2
C ... for snot
C     double precision idx(100),wk(300),r2(100),bes(0:20),han(0:20)

      if (el < -1d-12) call rx('us2hj not setup for linked basis')

      er2 = e*rmt*rmt
      call bessl(er2,lmax+1,phi,psi)
      rl = 1d0/rmt
      do  10  l = 0, lmax
        il = l
        if (job == 1) il = 0
        rl = rl*rmt
        vj = rl*phi(l)
        sj = rl*(l*phi(l)-er2*phi(l+1))/rmt
        vh = psi(l)/(rl*rmt)
        sh = (l*psi(l)-psi(l+1))/(rl*rmt*rmt)
        i1 = 1
        if (rofi(1) < 1d-10) i1 = 2
        do  20  ir = i1, nr
        b(ir,l)  = (vj*ul(ir,l) + sj*sl(ir,l))/rofi(ir)**(il+1)
   20   xi(ir,l) = (vh*ul(ir,l) + sh*sl(ir,l))/rofi(ir)**(il+1)
        if (i1 == 2) then
          xi(1,l) = 0
          xi(1,0) = xi(2,0)
        endif
   10 continue

      end
