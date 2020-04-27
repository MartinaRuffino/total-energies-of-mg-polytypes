      subroutine us2h(e,lmax,el,tl,rmt,rofi,ul,sl,nrx,nr,xi)
C- Linear combinations of (u,s) that match differentiably onto e
C  xi and b: true radial parts/r**l so the solid hankel functions are
C  hl(ilm) = xi(l)*cy(ilm)*yl(ilm) or jl(ilm) = b(l)*cy(ilm)*yl(ilm)
      implicit none
      integer nrx,nr,lmax
      double precision xi(nrx,0:1),ul(nr,0:1),sl(nr,0:1),
     .  rofi(1),el,tl,e,rmt
      integer ir,l,i1
      double precision phi(0:20),psi(0:20),sh,sj,vh,vj,rl,er2
C ... for snot
      double precision idx(100),wk(300),r2(100)

      if (el < -1d-12) call rx('us2h not setup for linked basis')

      er2=e*rmt*rmt
      call bessl(er2,lmax+1,phi,psi)
      rl = 1d0/rmt
      do  10  l = 0, lmax
        rl = rl*rmt
        vj = rl*phi(l)
        sj = rl*(l*phi(l)-er2*phi(l+1))/rmt
        vh = psi(l)/(rl*rmt)
        sh = (l*psi(l)-psi(l+1))/(rl*rmt*rmt)

C        print 333, vj,sj,vh,sh
C  333   format(4f12.6)

        i1 = 1
        if (rofi(1) < 1d-10) i1 = 2
        do  20  ir = i1, nr
   20   xi(ir,l) = (vh*ul(ir,l) + sh*sl(ir,l))/rofi(ir)**(l+1)
        if (i1 == 2) then
          xi(1,l) = 0
          xi(1,0) = xi(2,0)
        endif
   10 continue


C      do  50  ir = 1, nr
C   50 r2(ir) = rofi(ir)**2
C If you uncomment this check the call to hansr, it's changed!
C      call hansr(.5d0,0,lmax,1,lmax,e,r2,nrx,nr,idx,wk,0,xi)
C      call hansr(.3d0,0,lmax,1,lmax,e,r2,nrx,nr,idx,wk,0,xi)

C      print *, 'lmax,r2(nr),xi=',nr,r2(nr),xi(nr,2)
C      pause 'us2h'


      end
