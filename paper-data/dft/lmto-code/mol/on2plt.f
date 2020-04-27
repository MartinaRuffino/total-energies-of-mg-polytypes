      subroutine on2plt(ep1,ep2,lp1,lp2,nx,lmax,cof,
     .  rmt,xie,phi1,phi2,r,w,nr,n1,errs)
C- made from on2err, plot fit and product
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(nr),w(nr),fac1(0:10),fac2(0:10),
     .  cof(0:lmax,nx,0:lp1,0:lp2),errs(1),
     .  phi1(nr,0:1),phi2(nr,0:1),xie(nr,0:lmax,1),
     .  rout(901),pout(901),fout(901,5)
      integer fopna

      call getpr(ipr)
          rtop=2.5d0
      do 22 ir=1,n1
  22  if(r(ir) <= rtop) nrx=ir
      n2=nr-n1
      write(6,*) 'on2plt: enter l1,l2'
      read (5,*)  l01,l02
      ifi = fopna('pl',49,0)
      rewind 49
      fpi=16.d0*datan(1.d0)
      lsym=0
      if(lp1 == lp2.and.dabs(ep1-ep2) < 1d-10) lsym=1
c ------ get avg values on sphere to scale errors -------
      call hrmsav(rmt,lp1,ep1,fac1)
      call hrmsav(rmt,lp2,ep2,fac2)
      if(ipr >= 40) write(6,330) 1,(fac1(l),l=0,lp1)
      if(ipr >= 40) write(6,330) 2,(fac2(l),l=0,lp2)
  330 format(' rms avg of phi',i1,' on sphere: ',6f9.4)
c ------ start loop over pairs phi1,phi2 ------
      nout=0
      do 30 l1=0,lp1
      ltop=lp2
      if(lsym == 1) ltop=l1
      do 30 l2=0,ltop
      lm1=iabs(l1-l2)
      lm2=l1+l2
      do 30 lm=lm1,lm2,2
      if(l1 /= l01.or.l2 /= l02) goto 30
c ------ loop over points, make product and fit -----------
      fac=1d0/(fac1(l1)*fac2(l2)*fpi)
      nout=nout+1
      do 35 ir=1,nr
      rr=r(ir)
      prd=phi1(ir,l1)*phi2(ir,l2)*fac
      fit=0.d0
      do 10 ix=1,nx
   10 fit = fit + cof(lm,ix,l1,l2)*xie(ir,lm,ix)*fac
      jr=0
      if(ir <= nrx) jr=n2+ir
      if(ir > n1)  jr=ir-n1
c     print *, ir, nout, jr, r(ir), prd, fit
      if(jr > 0) then
        rout(jr)=rr
        pout(jr)=prd
        fout(jr,nout)=fit
        endif
  35  continue
  30  continue
c ------ write out --------------
      do 40 jr=1,n2+nrx
C|    write( 6,721) rout(jr),pout(jr),(fout(jr,j),j=1,nout)
C   40 write(49,721) rout(jr),pout(jr),(fout(jr,j),j=1,nout)
   40 write(49,721) rout(jr),
     .    min(pout(jr),5d0),(min(fout(jr,j),5d0),j=1,nout)
  721 format(6f12.6)

      stop
      end
