      subroutine ftcgen(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,d,nb1,nb2,xp,zp,wp,np,npgo,
     .  zc1,zc2,ri1,ri2,ifi,lchk,err,lpr)
C- Makes 2-center fit for distance d and one pair of phi-energies
C  lscal=1: improves condition of overlap matrix by scaling xi's
      implicit real*8 (a-h,p-z), integer (o)
      dimension ex1(8),ex2(8),err(100),cx(300),xp(np),zp(np),wp(np),
     .   ndim(0:20),nrhs(0:20),irhs(0:20),lx1(8),lx2(8),lx(20)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call sxlmnc(cx,10)
      lscal=1
      if(iprint() >= 30) write(6,550) d,np,npgo,lscal
  550 format(/' ftcgen:  d=',f10.6,'    np,npgo=',2i7,'    lscal=',i1)

C ------ setup ---
      call ftcdim(lx1,ex1,nx1,lx2,ex2,nx2,mmax,lmax,ep1,lp1,
     .   ep2,lp2,rmt1,rmt2,rsm1,rsm2,rsmp1,rsmp2,lsym,ndim,nrhs,
     .   ndimx,nsdmx,nrhsx,ncof,nxi,000)
      call defrr(ob,     ndimx*nrhsx*(mmax+1))
      call defrr(os,     nsdmx*(mmax+1))
      call dpzero(w(ob), ndimx*nrhsx*(mmax+1))
      call dpzero(w(os), nsdmx*(mmax+1))

C ------ allocate arrays for xi's and phi's -------
      npdim=min0(np,npgo)
      nlm = ((lmax+1)*(lmax+2))/2
      nlm1 = ((lp1+1)*(lp1+2))/2
      nlm2 = ((lp2+1)*(lp2+2))/2
      call defrr(oxi,    npdim*nlm*(nx1+nx2))
      call defrr(oph1,   npdim*nlm1)
      call defrr(oph2,   npdim*nlm2)
      call defrr(owk,    npdim)

c ------ loop over sets of npgo points; accumulate s,b ------
      ngo=(np-1)/npgo+1
CL      if(lpr == 1) write(71,760) lscal,np,npgo,ngo
  760 format(' lsc',i4,'   np',i7,'   npgo',i7,'   passes:',i4)
      do 40 igo=1,ngo
      ip=(igo-1)*npgo+1
      mp=npgo
      if(igo == ngo) mp=np-ip+1
      if(ipr >= 41) write(6,851) igo,mp,ip,ip+mp-1
  851 format(10x,'pass',i3,' :',i7,'  points,',i6,'   to',i6)
c ... set up phi's and xi's
      call ftcxip(lx1,ex1,nx1,lx2,ex2,nx2,rmt1,rmt2,rsmp1,rsmp2,
     .  rsm1,rsm2,lp1,ep1,lp2,ep2,cx,d,lscal,lmax,nx,lx,zc1,zc2,
     .  mp,xp(ip),zp(ip),wp(ip),w(oph1),w(oph2),w(oxi))
c ... add to b and s
      call ftcnrm(mmax,lmax,lx,nx,w(owk),wp(ip),mp,w(oxi),nlm,
     .  cx,w(oph1),w(oph2),lp1,lp2,lsym,w(os),w(ob),
     .  ndimx,nsdmx,nrhsx)
  40  continue
      call rlse(oxi)

C ------ solve the least-squares problem -----
      call ftcsev(mmax,ndim,nsdmx,w(os),lpr)
      call ftcslv(mmax,w(os),w(ob),ndim,nrhs,ndimx,nsdmx,nrhsx)
      if(lscal == 1)
     . call ftcusc(mmax,ndim,nrhs,nx1,lx1,ex1,nx2,lx2,ex2,
     .    rmt1,rmt2,w(ob),ndimx,nrhsx)
      call rlse(os)

c ------ sort coefficients into cof, write to file ifi ----
      call defrr(ocof,    ncof)
      call tcfsrt(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,ndim,irhs,
     .   w(ob),ndimx,nrhsx,w(ocof))
      if(ifi > 0) call dpdump(w(ocof),ncof,-ifi)

C|       call wreit(w(ocof),ncof,d)

C ------ check fit (use coarser mesh) ------
      if(lchk == 1) then
      mb1=15
      mb2=21
      call ftcchk(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,rmt1,rmt2,
     .  rsmp1,rsmp2,rsm1,rsm2,d,mb1,mb2,ri1,ri2,w(ocof),err)
      endif

      call rlse(ob)
      end
        subroutine wreit(cof,ncof,d)
        implicit real*8 (a-h,p-z)
        dimension cof(ncof)
        write(15,100) ncof,d
  100   format(i10,f10.5)
        write(15,101) cof
  101   format(1p,5d12.5)
        end
