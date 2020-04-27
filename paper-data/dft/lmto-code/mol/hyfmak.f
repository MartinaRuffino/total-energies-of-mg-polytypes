      subroutine hyfmak(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,np1,lp2,ep2,np2,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,nr1,nth1,w1,nr2,nth2,w2,
     .  d0,ndist,adec,a1,nalf,nb1,nb2,ri1,ri2,ifi)
C  Makes global 2-center expansion by fitting PHI1*PHI2 as
C  linear combination of XI's for NDIST distances simultaneously.
C  The PHI's and XIs are smoothed Hankel functions,
C  phi's given by RSMP1,RSMP2 and xi's by RSM1,RSM2
C  LP1,EP1 and LP2,EP2 describe the orbitals to be expanded.
C  LX1,EX1,NX1 and LX2,EX2,NX2 describe the XI-basis.
C  RI1,RI2,Nb1,Nb2,NPMX determine the bispherical integration mesh,
C  so the fit is done for the volume outside RI1 and RI2.
C  Rmt1,Rmt2 are the muffin-tin radii for which the fit is designed.
C  They are carried along to help select the correct table later on,
C  but they have no influence on the fit itself. They do enter into
C  the expressions for the fit error.
C  The fit distances are log distributed with ADEC stating at D0.
C  Fit coefficients have distance dependence exp(-M*A1), A1=1..NALF.
C  Output: TCF info and data on file IFI.
C  Bugs: hyfdim's check for symmetric case is incomplete.
      implicit none
      integer ndim(0:20),nrhs(0:20),ndimb(0:20),lx1(8),lx2(8),lx(20)
      integer lp1,lp2,lsym,mmax,nb1,nb2,nalf,ndimbx,ndimx,ndist,nrhsx,
     .  nr1,nr2,nsdmbx,nsdmx,nth1,nth2,nx1,nx2,np1,np2,ifi
      double precision ep1,ep2,rmt1,rmt2,zc1,zc2,ri1,ri2,
     .  rsm1,rsm2,rsmp1,rsmp2,w1,w2,
     .  ex1(8),ex2(8),dist(80),wist(80),err(20),cx(300),
     .  a1,adec,d0
      integer oph1,oph2,ob,of,os,ocof,obot,otop,oba,obb,
     .  oerr1,oerr2,oprd,osb,oxi,owp,oxp,ozp,i,idist,
     .  ipr,iprint,irep,lmax,m,nbisi,ncof,nlm,np,
     .  nrep,nx,nxi
      double precision d
      real w(1)
      common /w/ w

C --- Setup ---
      call getpr(ipr)
      call sxlmnc(cx,10)
      if (ndist > 80) call rx('HYFMAK: ndist too big')
      nbisi = 2*nb1*nb2
      np = nbisi + nr1*nth1 + nr2*nth2
      mmax = lp1+lp2
      stop 'update call to hyfdm'
C      call hyfdm(mmax,lp1,ep1,np1,lx1,ex1,nx1,lp2,ep2,np2,lx2,ex2,nx2,
C     .   rmt1,rmt2,rsm1,rsm2,lsym,ndim,nrhs,ndimx,nsdmx,nrhsx,ncof,
C     .   nxi,nalf,ndimb,ndimbx,nsdmbx)
      call defrr(of,   ndist*nalf)
      call hyfstf(adec,a1,d0,dist,wist,w(of),ndist,nalf)
C --- Write to log ---
      if (ipr >= 20) write(6,408) nb1,nb2,nr1,nth1,nr2,nth2,np,ri1,ri2
  408 format(' hyfmak: nbisl=',2i3,'   nr1,nth1=',2i3,
     .  '   nr2,nth2=',2i3,'    np=',i4,'    ri=',2f7.3)
CL      write(71,710) d0,ndist,rmt1,rmt2,adec,a1,nalf
  710 format(' ------- HYFMAK:  d0',f9.5,'   ndist',i3,'   r',2f9.4
     .   /' adec',f7.3,'   a1',f7.3,'   nalf',i3)
CL      write(71,991) mmax,(ndim(m),m=0,mmax)
  991 format(' mmax',i3,'   ndim',10i4)
CL      write(71,993) nb1,nb2,np,ri1,ri2
  993 format(' nb1,nb2',2i4,'   np',i5,'   ri1,ri2',2f7.3)
CL      write(71,990) 1,(lx1(i),i=1,nx1),(200,i=nx1+1,6),
CL     .   (ex1(i),i=1,nx1)
CL      write(71,990) 2,(lx2(i),i=1,nx2),(200,i=nx2+1,6),
CL     .   (ex2(i),i=1,nx1)
CL      write(71,992) lp1,lp2,ep1,ep2,rsm1,rsm2,lsym
  990 format(' xi',i1,'   l= ',6(1x,i1),'   e=',6f7.2)
  992 format(' phi   lmx=',2i2,'   e=',2f7.2,'   rsm=',2f7.3,
     .   '   lsym=',i1)
C --- Make bb and ss ---
      call defrr(obb,    ndimbx*nrhsx*(mmax+1))
      call defrr(osb,    nsdmbx*(mmax+1))
      call defrr(ob,     ndimx*nrhsx*(mmax+1))
      call defrr(os,     nsdmx*(mmax+1))
      call hyfint(mmax,lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
     .  ri1,ri2,rsmp1,rsmp2,rsm1,rsm2,
     .  nr1,nth1,w1,nr2,nth2,w2,
     .  lp1,ep1,lp2,ep2,lsym,w(os),w(ob),ndim,nrhs,ndimx,nsdmx,
     .  nrhsx,w(osb),w(obb),ndimbx,nsdmbx,dist,wist,w(of),ndist,nalf)
      call rlse(ob)
C --- Solve big lsqfit problem ---
      if (iprint() >= 60) call tm('start tcfslv')
      call tcfslv(mmax,w(osb),w(obb),ndimb,nrhs,ndimbx,nsdmbx,nrhsx)
      call rlse(osb)
C --- Setup for checking fit ---
      if (iprint() >= 60) call tm('start tcfchk')
      nrep = 3
      if (ndist == 1) nrep = 1
      do  70  irep = 1, nrep
      if (irep == 1) idist = 1
      if (irep == 2) idist = (ndist+1)/2
      if (irep == 3) idist = ndist
      d = dist(idist)
      call defrr(ob,     ndimx*nrhsx*(mmax+1))
      call hyfip1(d,d0,a1,mmax,ndim,nrhs,nalf,w(obb),w(ob),
     .  ndimx,ndimbx,nrhsx)
      stop 'update call to hyfxim'
C      call hyfxim(lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,
C     .  ri1,ri2,rsmp1,rsmp2,rsm1,rsm2,
C     .  nr1,nth1,w1,nr2,nth2,w2,
C     .  np,lp1,ep1,lp2,ep2,
C     .  cx,d,0,
C     .  lmax,nx,lx,zc1,zc2,oph1,oph2,owp,oxi,oxp,ozp)
      nlm = ((lmax+1)*(lmax+2))/2

C --- Check fit ---
      call defrr(oerr1,      nrhsx*(mmax+1))
      call defrr(oerr2,      nrhsx*(mmax+1))
      call defrr(oprd,       nrhsx*(mmax+1))
      call defrr(obot,       nrhsx*(mmax+1))
      call defrr(otop,       nrhsx*(mmax+1))
      call tcfchk(mmax,lx1,nx1,lx2,nx2,w(oxp),w(ozp),w(owp),np,
     .  nx,lx,w(oxi),nlm,lmax,w(oph1),w(oph2),ri1,ri2,w1,w2,
     .  rmt1,rmt2,zc1,zc2,lp1,ep1,lp2,ep2,lsym,w(ob),ndim,
     .   nrhs,ndimx,nrhsx,w(oerr1),w(oerr2),w(oprd),w(obot),w(otop),err)
      call rlse(ob)
   70 continue
      if (iprint() >= 60) call tm('done tcfchk')
C --- Sort coeffs into cof by ialf-columns ---
      call defrr(ocof,   ncof*nalf)
      call defrr(oba,    ndimx*nrhsx*(mmax+1)*nalf)
      call hyfsrt(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,ndim,nrhs,
     .  w(oba),w(obb),ndimx,nrhsx,ndimbx,nalf,w(ocof),ncof)
C --- Output into file ifi ---
      call hyfout(rmt1,rmt2,rsm1,rsm2,d0,adec,a1,ndist,nalf,mmax,ncof,
     .   nb1,nb2,ri1,ri2,nxi,lp1,ep1,lp2,ep2,lsym,lx1,ex1,nx1,
     .   lx2,ex2,nx2,w(ocof),err,ifi)
      call rlse(of)

      end
