      subroutine ftcmak(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d0,ndist,adec,nalf,
     .  nb1,nb2,ri1,ri2,jfi,ifi)
C- read pointwise 2c-fits, make chebyshev coeffs, output final fit
      implicit real*8 (a-h,p-z), integer (o)
      dimension lx1(8),lx2(8),lx(20),ndim(0:20),err(100),
     .  ex1(8),ex2(8),dist(120),nrhs(0:20),irhs(0:20)
      real w(1)
      common /w/ w
C --- setup ---
      call getpr(ipr)
      if (ndist > 120) call rx('ftcmak: ndist too big')
      call ftcdis(d0,ndist,adec,dist)
      call ftcdim(lx1,ex1,nx1,lx2,ex2,nx2,mmax,lmax,ep1,lp1,
     .  ep2,lp2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym,ndim,nrhs,
     .  ndimx,nsdmx,nrhsx,ncof,nxi,000)

C --- make chebyshev coeffs ---------
      if(ndist == 1.and.nalf /= 1) write(6,*) 'ftcmak:  nalf set to 1'
      if(ndist == 1) nalf=1
      call defrr(ocof,   ncof*ndist)
      call defrr(oche,   ncof*nalf)
      call ftcinc(ncof,ndist,w(ocof),jfi)
      call ftcmch(ncof,ndist,nalf,w(ocof),w(oche))

C --- assemble fit, compare to cof ---------
      call defrr(ocuf,   ncof*ndist)
      call ftcass(ncof,ndist,nalf,dist,adec,d0,w(oche),w(ocuf))
      call ftccmp(ncof,ndist,w(ocuf),w(ocof),erravg,errmax,j0)
      if(ipr >= 20) write(6,330) errmax,j0,erravg
  330 format(10x,'max cheby fit error=',f9.3,' %   for i=',i6/
     .       10x,'avg cheby fit error=',f9.3,' %')
CL      write(71,710) nalf,erravg,errmax,j0
  710 format(' nalf=',i3,'    cheb err:  avg',f7.2,' %   top',
     .  f7.2,' %  at',i5)
      call rlse(ocuf)

c --- output into file ifi -------------------------
      call dpzero(err,100)
      call hyfout(rmt1,rmt2,rsm1,rsm2,d0,adec,adec,ndist,nalf,mmax,
     .   ncof,nb1,nb2,ri1,ri2,nxi,lp1,ep1,lp2,ep2,lsym,lx1,ex1,nx1,
     .   lx2,ex2,nx2,w(oche),err,ifi)

c --- get error of final tables at d0 ---
      d=d0
      x=dexp(adec*(d0-d))
      y=2d0*x-1d0
      call defrr(owk,  2*ncof)
      call ropecs(y,nalf,ncof,w(owk),w(oche),w(ocof))
        call rlse(oche)
      mb1=15
      mb2=21
      call ftcchk(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,rmt1,rmt2,
     .  rsmp1,rsmp2,rsm1,rsm2,d,mb1,mb2,ri1,ri2,w(ocof),err)
      call rlse(ocof)

      end
