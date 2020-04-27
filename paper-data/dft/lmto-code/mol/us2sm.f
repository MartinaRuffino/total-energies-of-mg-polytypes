      subroutine us2sm(rmt,a,nr,rus,wus,ul,sl,lmxa,npwr,nrad,cof,
     .  rad,wrad,ul2,sl2)
C- Makes smoothed polynomical approximation to u,s
      implicit none
      integer nr,npwr,nrad,lmxa
      double precision rmt,a,rus(1),wus(1),rad(1),wrad(1),
     .  ul(nr,0:1),sl(nr,0:1),ul2(nrad,0:1),sl2(nrad,0:1)
      integer ip,jp,ir,is,nchd,l,nc,npmx,nnu,nns,nrmx
      parameter (npmx=15)
      double precision sumu,suml,eru,ers,
     .  cof(npwr,2,0:1),s(npmx**2),wk(npmx),wk2(3),q(2,2),p(npmx*2)

C --- Setup ----
      if (npwr > npmx) call rx('us2sm:  npwr gt npmx')
C     call junk2(rus,sl,nr,1)
      call mklegw(nrad,rad,wrad,0)
      do  5  ir = 1, nrad
      rad(ir) = rmt/2*rad(ir) + rmt/2
    5 wrad(ir) = wrad(ir)*rmt/2

C --- Fit ul,sl, one l at at time ---
      do  80  l = 0, lmxa

C --- Count number of nodes in ul,sl ---
      nnu = 0
      do  10  ir = 3, nr
   10 if (ul(ir,l)*ul(ir-1,l) < 0) nnu = nnu+1

C --- Normal matrix and rhs for LS fit ---
      is = 0
      do  20  ip = 1, npwr
        do  25  jp = 1, ip
          is = is+1
          s(is) = 0
          do  27  ir = 1, nr
   27     s(is) = s(is) + wus(ir)*rus(ir)**(ip+jp+2*l+2-2)
   25   continue

        sumu = 0d0
        do  28  ir = 1, nr
   28   sumu = sumu + wus(ir)*rus(ir)**(ip+l+2-1-1)*ul(ir,l)
        cof(ip,1,l) = sumu
        suml = 0d0
        do  29  ir = 1, nr
   29   suml = suml + wus(ir)*rus(ir)**(ip+l+2-1-1)*sl(ir,l)
        cof(ip,2,l) = suml
   20 continue
      call chlr2f(s,wk,npwr,nchd)
      if (nchd < npwr) call rx('US2SM: normal matrix not pos def')

C --- Fit with constraint for val, slo at rmt ---
      do  30  ip = 1, npwr
      p(ip+npwr) = (ip+l-1)*rmt**(ip+l-2)
   30 p(ip) = rmt**(ip+l-1)
      q(1,1) = 1
      q(2,1) = 0
      q(1,2) = 0
      q(2,2) = 1
      nc = 2
      call lsqfc1(s,p,npwr,nc,wk2,wk)
      call lsqfc2(cof(1,1,l),q,s,p,npwr,2,nc,wk2)

C --- Check errors ---
C      rewind 80
C      write(80,*) nr-1, 5
      eru = 0
      ers = 0
      do  40  ir = 2, nr
        sumu = 0
        do  42  ip = 1, npwr
   42   sumu = sumu + cof(ip,1,l)*rus(ir)**(ip+l-1)
        eru = eru + wus(ir)*rus(ir)**2*(sumu-ul(ir,l)/rus(ir))**2
        suml = 0
        do  43  ip = 1, npwr
   43   suml = suml + cof(ip,2,l)*rus(ir)**(ip+l-1)
        ers = ers + wus(ir)*rus(ir)**2*(suml-sl(ir,l)/rus(ir))**2
C        write(80,333) rus(ir), ul(ir,l)/rus(ir), sumu,
C     .    sl(ir,l)/rus(ir), suml
C  333   format(f8.5,4f12.6)
   40 continue

      print 345, l, nnu, dsqrt(eru/nr), dsqrt(ers/nr)
  345 format(' us2sm:  l=',i1,'  nnode=',i1,'  rms err (u,s)=',2f10.6)

   80 continue
      call us2sm2(npwr,nrad,lmxa,cof,rad,ul2,sl2)
      end
      subroutine us2sm2(npwr,nrad,lmxa,cof,rad,ul,sl)
C- Evaluates (u,s) on a mesh from a power series expansion
      implicit none
      integer npwr,nrad,lmxa
      double precision ul(nrad,0:1),sl(nrad,0:1),cof(npwr,2,0:1),rad(1)
      integer ip,ir,l
      double precision suml,sumu

      do  80  l = 0, lmxa
        do  50  ir = 1, nrad
          sumu = 0
          do  52  ip = 1, npwr
   52     sumu = sumu + cof(ip,1,l)*rad(ir)**(ip+l-1+1)
          ul(ir,l) = sumu
          suml = 0
          do  53  ip = 1, npwr
   53     suml = suml + cof(ip,2,l)*rad(ir)**(ip+l-1+1)
          sl(ir,l) = suml
   50   continue

   80 continue

C      do  60  ir = 1, nrad
C   60 write(*,333) rad(ir), (ul(ir,l)/rad(ir), l=0,lmxa),
C     .    (sl(ir,l)/rad(ir), l=0,lmxa)
C  333 format(f8.5,6f12.6)
C      pause

      end
