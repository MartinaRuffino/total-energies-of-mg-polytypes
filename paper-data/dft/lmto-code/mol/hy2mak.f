      subroutine hy2mak(lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,lmxl,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,
     .  nr1,nth1,rad1,wrad1,w1,phr1,nr2,nth2,rad2,wrad2,w2,phr2,
     .  d0,ndist,adec,a1,nalf,nb1,nb2,ifi)
C- Makes global 2-center expansion by fitting PHI1*PHI2 as a
C  linear combination of XI's for NDIST distances simultaneously.
C  The XI's are smoothed Hankel functions, described by NX,LX,EX,RSM
C  for sites 1 and 2.  Any EX>=0 is interpreted as a Gaussian with
C  width = EX (EX=0 makes gaussian with width=RSM).
C  The PHI's are the orbitals to be expanded.  Outside RMT they are
C  smoothed Hankels, described by RSMP,LP,EP for sites 1 and 2.
C  Inside RMT, the phi's PH are passed as a tabulated function on
C  radial mesh RM and weights W with NR points for sites 1 and 2.
C  Fit outside RMT is done on a bispherical integration mesh;
C  NB1,NB2 determine the mesh.  Fit inside RMT is done on the tabulated
C  radial mesh RM and weights W, expanded into NR*NTH points for 2d.
C  In total, 2*NB1*NB2 + NR1*NTH1 + NR2*NTH2 points are generated.
C  The fit distances are log distributed with ADEC stating at D0.
C  Fit coefficients have distance dependence exp(-M*A1), A1=1..NALF.
C  Output: TCF info and data on file IFI.
C  Bugs: hyfdim's check for symmetric case was based on an old hyfmak,
C  and is incomplete.
      implicit none
      integer ndim(0:20),nrhs(0:20),ndimb(0:20),lx1(8),lx2(8),lx(20)
      integer lp1,lp2,lsym,mmax,nb1,nb2,nalf,ndimbx,ndimx,ndist,nrhsx,
     .  nr1,nr2,nsdmbx,nsdmx,nth1,nth2,nx1,nx2,npmx,lmxl,ifi
      double precision ep1,ep2,rmt1,rmt2,zc1,zc2,
     .  rad1(1),rad2(1),wrad1(1),wrad2(1),rsm1,rsm2,rsmp1,rsmp2,w1,w2,
     .  ex1(8),ex2(8),dist(80),wist(80),err(20),cx(300),
     .  phr1(nr1,0:1),phr2(nr2,0:1),
     .  a1,adec,d0
      integer oph1,oph2,ob,of,os,ocof,obot,oidx,otop,oba,obb,
     .  oerr1,oerr2,osb,oxi,owk,oyl,owp,oxp,ozp,i,idist,ioff1,ioff2,
     .  ipr,iprint,irep,ix,lmax,m,nadd,nbisi,ncof,nlm,nlmp1,nlmp2,np,
     .  nrep,nx,nxi
      double precision d
      real w(1)
      common /w/ w

C --- Setup ---
      call getpr(ipr)
      call sxlmnc(cx,10)
      if (ndist > 80) call rx('HY2MAK: ndist too big')
      nbisi = 2*nb1*nb2
      npmx = nbisi + nr1*nth1 + nr2*nth2
      mmax = lp1+lp2
      call hyfdim(mmax,lx1,ex1,nx1,lx2,ex2,nx2,lp1,ep1,lp2,ep2,
     .   rmt1,rmt2,rsm1,rsm2,lsym,ndim,nrhs,ndimx,nsdmx,nrhsx,ncof,
     .   nxi,nalf,ndimb,ndimbx,nsdmbx)
      call defrr(of,   ndist*nalf)
      call hyfstf(adec,a1,d0,dist,wist,w(of),ndist,nalf)
      lmax = 0
      do  10  ix = 1, nx1
   10 lmax = max(lmax,lx1(ix))
      do  11  ix = 1, nx2
   11 lmax = max(lmax,lx2(ix))
      nlm = ((lmax+1)*(lmax+2))/2
C --- Write to log ---
      if (ipr >= 20) write(6,408) nb1,nb2,nr1,nth1,nr2,nth2,npmx
  408 format(' hy2mak: nbisi=',2i3,'   nr1,nth1=',2i3,
     .  '   nr2,nth2=',2i3,'    np=',i4)
CL      write(71,408) nb1,nb2,nr1,nth1,nr2,nth2,npmx
CL      WRITE(71,710) D0,NDIST,RMT1,RMT2,ADEC,A1,NALF
  710 FORMAT(' ------- HYFMAK:  D0',F9.5,'   NDIST',I3,'   R',2F9.4
     .   /' ADEC',F7.3,'   A1',F7.3,'   NALF',I3)
CL      WRITE(71,991) MMAX,(NDIM(M),M=0,MMAX)
  991 FORMAT(' MMAX',I3,'   NDIM',10I4)
CL      WRITE(71,993) NB1,NB2,NPMX
  993 FORMAT(' NB1,NB2',2I4,'   NP',I5)
CL      WRITE(71,990) 1,(LX1(I),I=1,NX1),(200,I=NX1+1,6),
CL     .   (EX1(I),I=1,NX1)
CL      WRITE(71,990) 2,(LX2(I),I=1,NX2),(200,I=NX2+1,6),
CL     .   (EX2(I),I=1,NX1)
CL      WRITE(71,992) LP1,LP2,EP1,EP2,RSM1,RSM2,LSYM
  990 FORMAT(' XI',I1,'   L= ',6(1X,I1),'   E=',6F7.2)
  992 FORMAT(' PHI   LMX=',2I2,'   E=',2F7.2,'   RSM=',2F7.3,
     .   '   LSYM=',I1)
C --- Make bb and ss ---
      call defrr(obb,    ndimbx*nrhsx*(mmax+1))
      call defrr(osb,    nsdmbx*(mmax+1))
      call defrr(ob,     ndimx*nrhsx*(mmax+1))
      call defrr(os,     nsdmx*(mmax+1))
      call hy2int(mmax,lx1,ex1,nx1,lx2,ex2,nx2,nb1,nb2,lmxl,
     .  rmt1,rmt2,rsmp1,rsmp2,rsm1,rsm2,
     .  nr1,nth1,rad1,wrad1,w1,phr1,nr2,nth2,rad2,wrad2,w2,phr2,
     .  npmx,lp1,ep1,lp2,ep2,lsym,w(os),w(ob),ndim,nrhs,ndimx,nsdmx,
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
      call defrr(oxp,    npmx)
      call defrr(owp,    npmx)
      call defrr(ozp,    npmx)
      call hy2msh(nr1,nth1,rad1,wrad1,w1,nr2,nth2,rad2,wrad2,w2,nadd,
     .  rmt1,rmt2,d,nb1,nb2,w(oxp),w(ozp),w(owp),npmx,nbisi,np,zc1,zc2)
      if (np /= npmx) call rx('hy2mak: mismatch npmx, np')
      call defrr(oxi,    np*nlm*(nx1+nx2))
      nlmp1 = ((lp1+1)*(lp1+2))/2
      call defrr(oph1,   np*nlmp1)
      nlmp2 = ((lp2+1)*(lp2+2))/2
      call defrr(oph2,   np*nlmp2)
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm,nlmp1,nlmp2))
      ioff1 = nbisi
      ioff2 = nbisi + nr1*nth1

C --- Tabulate xi's and phi's for this mesh ----------
      call defrr(owk,    np)
      call defrr(oidx,   np*2)
      call defrr(oyl,    np*max(6,nlm))
      stop 'hy2mak: update call to hy2xi'
C      call hy2xi(np,w(oxp),w(ozp),w(owp),0,lmax,nlm,
C     . zc1,nx1,lx1,ex1,nlmp1,rsmp1,rsm1,lp1,ep1,phr1,ioff1,nr1,nth1,
C     . zc2,nx2,lx2,ex2,nlmp2,rsmp2,rsm2,lp2,ep2,phr2,ioff2,nr2,nth2,
C     . cx,w(oidx),w(owk),w(oyl),w(oph1),w(oph2),nx,lx,w(oxi))
      call rlse(owk)
C --- Check fit ---
      call defrr(oerr1,      nrhsx*(mmax+1))
      call defrr(oerr2,      nrhsx*(mmax+1))
      call defrr(obot,       nrhsx*(mmax+1))
      call defrr(otop,       nrhsx*(mmax+1))
      stop 'hy2mak: update call to tcfchk'
C      call tcfchk(mmax,lx1,nx1,lx2,nx2,w(oxp),w(ozp),w(owp),np,
C     .  nx,lx,w(oxi),nlm,lmax,w(oph1),w(oph2),rmt1,rmt2,w1,w2,
C     .  rmt1,rmt2,zc1,zc2,lp1,ep1,lp2,ep2,lsym,w(ob),ndim,
C     .   nrhs,ndimx,nrhsx,w(oerr1),w(oerr2),w(obot),w(otop),err)
      call rlse(ob)
   70 continue
      if (iprint() >= 60) call tm('done tcfchk')
C --- Sort coeffs into cof by ialf-columns ---
      call defrr(ocof,   ncof*nalf)
      call defrr(oba,    ndimx*nrhsx*(mmax+1)*nalf)
      call hyfsrt(mmax,lx1,nx1,lx2,nx2,lp1,lp2,lsym,ndim,nrhs,
     .   w(oba),w(obb),ndimx,nrhsx,ndimbx,nalf,w(ocof),ncof)
C --- Output into file ifi ---
      call hyfout(rmt1,rmt2,rsm1,rsm2,d0,adec,a1,ndist,nalf,mmax,ncof,
     .   nb1,nb2,rmt1,rmt2,nxi,lp1,ep1,lp2,ep2,lsym,lx1,ex1,nx1,
     .   lx2,ex2,nx2,w(ocof),err,ifi)
      call rlse(of)

      end
