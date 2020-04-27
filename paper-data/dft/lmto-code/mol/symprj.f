      subroutine symprj(ia,nlm,pos,nrc,g,ag,ng,cy,sym)
c  Makes symmetry projection matrices for one class (see iii,p40).
c  Matrices are those needed to make sym fct on atom ia.
      implicit real*8 (a-h,p-z), integer(o)
      dimension g(9,ng),ag(3,1),pos(3,nrc),cy(1),sym(nlm,nlm,nrc)
      real w(1)
      common /w/ w
      call getpr(ipr)
      lmax=ll(nlm)
      if(ipr >= 60) write(6,340) nrc,ia,lmax,ng
  340 format(/' make projectors for  nrc, ia, lmax, ng=',4i5)
      if(ipr >= 60) write(6,341) (ja,(pos(m,ja),m=1,3),ja=1,nrc)
  341 format('    site',i3,'    pos=',3f11.5)
      call dpzero(sym,    nlm*nlm*nrc)
      call defrr(oxmat,   nlm*nlm)
      call defrr(ormat,   nlm*nlm)
      call defrr(op,      3*nlm)
      call ylmrt0(lmax,nlm,w(oxmat),w(op),cy)
      wgt=1.d0/ng
      do 10 ig=1,ng
      call gpfnd2(g(1,ig),ag(1,ig),pos(1,ia),nrc,pos,ja)
      if(ja == 0) call rx('symprj: no site found mapped onto ia')
      call ylmrt1(lmax,nlm,g(1,ig),w(ormat),w(oxmat),w(op),cy)
  10  call dpadd(sym(1,1,ja),w(ormat),1,nlm*nlm,wgt)
      call rlse(oxmat)
c ------ this part to print out matrices -----------
      if(ipr < 80) return
      do 20 ja=1,nrc
      write(6,727) ia,ja
  727 format(/' symmetry matrices for   ia=',i3,'   ja=',i3)
      ilm=0
      do 21 l=0,lmax
      jlm1=l*l+1
      jlm2=(l+1)**2
      do 21 m=1,2*l+1
      ilm=ilm+1
  21  write(6,210) (sym(ilm,jlm,ja),jlm=jlm1,jlm2)
  210 format(1x,9f8.4)
  20  continue
      end
c ------- rotpnt -----------
      subroutine rotpnt(p,q,g)
      implicit real*8 (a-h,p-z), integer(o)
      dimension p(3),q(3),g(3,3),h(3)
      do 1 i=1,3
      h(i)=0.d0
      do 1 j=1,3
  1   h(i)=h(i)+g(i,j)*p(j)
      do 2 i=1,3
  2   q(i)=h(i)
      return
      end
