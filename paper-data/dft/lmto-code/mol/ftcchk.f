      subroutine ftcchk(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d,nb1,nb2,ri1,ri2,cof,err)
C- Checks 2cf at one distance d. Input is cof.
      implicit real*8 (a-h,p-z), integer (o)
      dimension ex1(8),ex2(8),err(100),cx(300),cof(1),
     .   ndim(0:20),nrhs(0:20),irhs(0:20),lx1(8),lx2(8),lx(20)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call sxlmnc(cx,10)
      if(ipr >= 30) write(6,300) nb1,nb2,d
  300 format(/' ftcchk:  nb1,nb2=',2i6,'   d=',f10.5)

C ------ setup ---
      call ftcdim(lx1,ex1,nx1,lx2,ex2,nx2,mmax,lmax,ep1,lp1,
     .   ep2,lp2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym,ndim,nrhs,
     .   ndimx,nsdmx,nrhsx,ncof,nxi,000)

C ------ make b from cof -------
      call defrr(ob,   ndimx*nrhsx*(mmax+1))
      call tcfusr(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,ndim,irhs,
     .   w(ob),ndimx,nrhsx,cof)

C ------ make mesh ----------
      npmx=2*nb1*nb2   +1000
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
C|    call bisinl(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
C|   .  zc1,zc2)
      call holint(ri1,ri2,d,nb1,nb2,w(oxp),w(ozp),w(owp),np,npmx,
     .  zc1,zc2)

C ------ allocate arrays for xi's and phi's -------
      nlm = ((lmax+1)*(lmax+2))/2
      nlm1 = ((lp1+1)*(lp1+2))/2
      nlm2 = ((lp2+1)*(lp2+2))/2
      call defrr(oxi,    np*nlm*(nx1+nx2))
      call defrr(oph1,   np*nlm1)
      call defrr(oph2,   np*nlm2)

C ------ check fit ---
      lscal=0
      call ftcxip(lx1,ex1,nx1,lx2,ex2,nx2,rmt1,rmt2,rsmp1,rsmp2,
     .  rsm1,rsm2,lp1,ep1,lp2,ep2,cx,d,lscal,lmax,nx,lx,zc1,zc2,
     .  np,w(oxp),w(ozp),w(owp),w(oph1),w(oph2),w(oxi))
      call defrr(oerr1,      nrhsx*(mmax+1))
      call defrr(oerr2,      nrhsx*(mmax+1))
      call defrr(owk,        np)
      call defrr(ofit,       np)
      call defrr(owerr,      np)
      call ftcerr(d,mmax,w(owp),np,w(oxp),w(ozp),zc1,zc2,w(owerr),
     .  w(owk),w(ofit),nx,lx,w(oxi),nlm,lmax,w(oph1),w(oph2),
     .  lp1,ep1,lp2,ep2,rmt1,rmt2,lsym,w(ob),ndimx,nrhsx,
     .  w(oerr1),w(oerr2),err)
      call rlse(ob)
      end
