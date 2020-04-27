      subroutine hy2msh(nr1,nth1,rad1,wrad1,w1,nr2,nth2,rad2,wrad2,w2,
     .  nadd,rmt1,rmt2,d,n1,n2,xp,zp,wp,npmx,nbisi,np,zc1,zc2)
C- Bisinl for plane integral; add points in plane inside spheres
      implicit none
      integer nr1,nth1,nr2,nth2,nadd,n1,n2,np,npmx,nbisi
      double precision rmt1,rmt2,d,xp(1),zp(1),wp(1),
     .  rad1(1),wrad1(1),rad2(1),wrad2(1),zc1,zc2,w1,w2
      integer iprint

      call pshpr(iprint()-30)
      call bisinl(rmt1,rmt2,d,n1,n2,xp,zp,wp,np,npmx,zc1,zc2)
      nbisi = np

      nadd = 0
      call hymsh2(nr1,nth1,rad1,wrad1,zc1,w1,np,nadd,xp,zp,wp)
      call hymsh2(nr2,nth2,rad2,wrad2,zc2,w2,np,nadd,xp,zp,wp)
      call poppr
      np = np+nadd
      if (np > npmx) call rx('hy2msh: np gt npmx')

      end
      subroutine hymsh2(nr,nth,rad,wrad,zc,w,np,nadd,xp,zp,wp)
C- Add points inside a sphere for 2d mesh
      implicit none
      integer nr,nth,np,nadd
      double precision zc,xp(1),zp(1),wp(1),rad(1),wrad(1),w
      double precision pi2
      integer nmx
      parameter (nmx=40)
      double precision cth(nmx),wcth(nmx)
      integer i,j,iprint

      if (nth > nmx) call rx('hy2sp1:  nth gt nmx')
      pi2 = 2*4*datan(1d0)
      call mklegw(nth,cth,wcth,iprint())
      do  1  j = 1, nth
      do  1  i = 1, nr
        nadd = nadd+1
        xp(np+nadd) = rad(i)*dsqrt(1-cth(j)**2)
        zp(np+nadd) = rad(i)*cth(j) + zc
        wp(np+nadd) = wrad(i)*w*wcth(j)*rad(i)**2*pi2
    1 continue

      end
