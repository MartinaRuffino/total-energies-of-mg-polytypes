      subroutine on2mak(llink,lx,ex,nx,lp1,ep1,add1,el1,
     .  lp2,ep2,add2,el2,rmt,rsm,rsmp,ri1,ri2,ifi)
C- One-center fit of smoothed, linked hankel products.  Fit for r>rmt.
C ----------------------------------------------------------------------
Ci Inputs
Ci   llink
Ci   lx    :product basis fitting l-cutoffs
Ci   ex    :product basis fitting energies
Ci   nx    :number of product basis fitting energies
Ci         :the following three entries correspond to
Ci   lp1   :lmto basis l-cutoff, first function
Ci   ep1   :lmto basis energy, first function
Ci   add1  :coefficient to linked part of basis, 1st function
Ci   el1   :?
Ci   lp2   :lmto basis l-cutoff, second function
Ci   ep2   :lmto basis energy, second function
Ci   add2  :coefficient to linked part of basis, 2nd function
Ci   el2   :?
Ci   rmt   :augmentation radius, in a.u., by species
Ci   rsm   :smoothing radius
Ci   rsmp  :smoothing radius for basis function
Ci   ri1   :not used
Ci   ri2   :not used
Ci   ifi   :file logical unit (for saving fit)
Co Outputs
Co   Fit to product basis is written to disk.
Cl Local variables
Cl   lmax  :maximum l in product basis
Cr Remarks
Cr
Cu Updates
Cu   24 Jun 04
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer llink,lp1,lp2,nx,ifi,lx(50)
      double precision add1(1),add2(1),el1,el2,ep1,ep2,ri1,ri2,ex(50),
     .  rmt,rsm,rsmp
C ... Local parameters
      integer nmax,ndef,stdo,jp(50),lxj(50),ipr,i,ir,j,lmax,lsym,lxi,
     .  nalf,ncof,ndim,nr,nglob
      parameter (nmax=401)
      double precision rofi(nmax),wofi(nmax),rsq(nmax),aa,tol,rmax,
     .  err(50)
      integer ob,ocof,oidx,ophi1,ophi2,ophil,os,owk,oxi,oxie
      integer w(1)
      common /w/ w

      stdo = nglob('stdo')
      call getpr(ipr)
      lmax = 0
      do  2  j = 1, nx
      lxj(j) = min0(lx(j),lp1+lp2)
    2 lmax = max0(lmax,lxj(j))
      ncof = (lmax+1)*nx*(lp1+1)*(lp2+1)
      call defrr(ocof,       -ncof)

C --- Setup for radial integration ---
      rmax = 30d0
      aa = 0.03d0
      call radint(rmt,aa,rmax,nmax,rofi,wofi,nr)
      do  54  ir = 1, nr
   54 rsq(ir) = rofi(ir)**2
      if (ipr >= 50) write(stdo,787) nr
  787 format(' on2mak: n1=',i5)

C --- Write to log ---
CL      write(71,710) rmt,rsm,rsmp,lmax
  710 format(' ------- on2mak:  r',f8.4,'  rsm',f8.4,
     .  '  rsmp',f8.4,'  lmax',i3)
CL      write(71,993) ri1,ri2,rmax,aa,nr
CL      write(71,996) (lx(i),i=1,nx),(200,i=nx+1,6),(ex(i),i=1,nx)
CL      write(71,997) lp1,lp2,ep1,ep2
  993 format(' ri',2f8.4,'    rmax',f8.2,'    a',f7.3,'   nr',i5)
  996 format(' xi:   l= ',6(1x,i1),'   e=',6f7.2)
  997 format(' phi   lmx=',2i2,'   e=',2f7.2)

C --- Radial parts of the xi and gaussians ---
      call defrr(oxie, -nr*(lmax+1)*nx)
      call defrr(oidx,  nr*2)
      call defrr(owk,   nr*6)
      tol = 1d-10
      call hansrg(rsm,lmax,nx,lxj,ex,rsq,nr,nr,tol,w(oidx),w(owk),001,
     .  w(oxie))
      call rlse(oidx)

C --- Radial parts of the phi1,phi2 ---
      call defrr(ophi1, nr*(lp1+1))
      call defrr(ophi2, nr*(lp2+1))
      call defrr(ophil, nr*max(lp1+1,lp2+1))
      call defrr(oidx,  nr*2)
      call defrr(owk,   nr*6)
      call on2phl(ep1,lp1,llink,rsmp,el1,add1,rsq,nr,
     .  w(oidx),w(owk),w(ophi1),w(ophil))
      call on2phl(ep2,lp2,llink,rsmp,el2,add2,rsq,nr,
     .  w(oidx),w(owk),w(ophi2),w(ophil))
      call rlse(ophil)

C --- Big loop over lxi ---
      do  80  lxi = 0, lmax
      ndim = 0
      do  1  j = 1, nx
      if (lxj(j) >= lxi) ndim = ndim+1
    1 if (lxj(j) >= lxi) jp(ndim) = j
      if (ndim == 0) goto 80
      call defrr(ob,     ndim*(lp1+1)*(lp2+1))
      call defrr(os,     ndim*ndim)
C --- Radial integrals for overlap and rhs ---
      call defrr(oxi,    nr*ndim)
      call on2int(rofi,wofi,nr,ndim,jp,lxi,w(oxie),lmax,w(oxi),
     .   lp1,lp2,w(ophi1),w(ophi2),w(os),w(ob))
      call rlse(oxi)
C --- Solve the least-squares problem ---
      call defrr(owk,    ndim)
      call chlr2f(w(os),w(owk),ndim,ndef)
      if (ndef /= ndim) call rx('on2mak:  matrix not pos definite''')
      call chlr2s(w(os),w(ob),ndim,(lp1+1)*(lp2+1))
C --- Copy to cof ---
      call hk0cop(ndim,jp,lxi,lp1,lp2,lmax,nx,w(ob),w(ocof))
      call rlse(ob)
   80 continue
C --- Fit errors ---
      call on2err(ep1,ep2,lp1,lp2,nx,lmax,-1,w(ocof),rmt,
     .  w(oxie),w(ophi1),w(ophi2),rofi,wofi,nr,err)
C     call on2plt(ep1,ep2,lp1,lp2,nx,lmax,w(ocof),rmt,
C    .  w(oxie),w(ophi1),w(ophi2),rofi,wofi,nr,nr,err)
C --- Output on file ifi ---
      i = 0
      nalf = 1
      lsym = 3
      call hyfout(rmt,rmt,rsm,rsm,0d0,0d0,0d0,i,nalf,lmax,
     .  ncof,nr,i,ri1,ri2,i,lp1,ep1,lp2,ep2,lsym,lxj,ex,nx,
     .  lxj,ex,nx,w(ocof),err,ifi)
      call rlse(ocof)

      end
      subroutine on2int(r,w,nr,ndim,jp,lxi,xie,lmax,xi,lp1,lp2,
     .  phi1,phi2,s,b)
C- (kernel called by on2mak) radial integrals for overlap and rhs
C ----------------------------------------------------------------------
Ci Inputs
Ci   r     :mesh of radial points
Ci   w     :integration weights for radial mesh
Ci   nr    :number of radial mesh points
Ci   ndim  :number of fitting functions for which to make integrals
Ci   jp    :for ith fitting function, sue xie(:,:,jp(i))
Ci         :(jp is permutation table)
Ci   lxi   :l cutoffs for each function in c.d. basis, by species
Ci   xie   :fitting functions from which list is taken; see jp
Ci   lmax  :maximum l for a given site
Ci   xi    :work array holding list of product functions
Ci   lp1   :form products phi1(1:lp1)*phi2(1:lp2) for pairs 1:lp1,1:lp2
Ci   lp2   :form products phi1(1:lp1)*phi2(1:lp2) for pairs 1:lp1,1:lp2
Ci   phi1  :form products phi1(1:lp1)*phi2(1:lp2) for pairs 1:lp1,1:lp2
Ci   phi2  :form products phi1(1:lp1)*phi2(1:lp2) for pairs 1:lp1,1:lp2
Co Outputs
Co   s     :overlap <xi_i xi_j>
Co   b     :integrals for rhs <phi1 phi2 xi_i>
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Jun 04
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,ndim,jp(1),lxi,lmax,lp1,lp2
      double precision r(nr),w(nr),xie(nr,0:lmax,1),
     .  phi1(nr,1),phi2(nr,2),s(1),b(ndim,lp1+1,lp2+1),xi(nr,1)
C ... Local parameters
      integer i,ir,is,ix,j,k1,k2
      double precision sum

C --- Assemble list of product functions from xie ---
      do  10  i = 1, ndim
      ix = jp(i)
      do  10  ir = 1, nr
   10 xi(ir,i) = xie(ir,lxi,ix)

C --- Overlap matrix ---
      is = 0
      do  20  i = 1, ndim
      do  20  j = 1, i
      is = is+1
      sum = 0
      do  21  ir = 1, nr
   21 sum = sum + xi(ir,i)*xi(ir,j)*w(ir)*r(ir)**2
   20 s(is) = sum

C --- Integrals for rhs ---
      do  30  i = 1, ndim
      do  30  k1 = 1, lp1+1
      do  30  k2 = 1, lp2+1
      sum = 0.d0
      do  31  ir = 1, nr
   31 sum = sum + phi1(ir,k1)*phi2(ir,k2)*xi(ir,i)*w(ir)*r(ir)**2
   30 b(i,k1,k2) = sum
      end
      subroutine hk0cop(ndim,jp,lxi,lp1,lp2,lmax,nx,b,cof)
      implicit none
C ... Passed parameters
      integer jp(1)
      integer ndim,lxi,lp1,lp2,lmax,nx
      double precision b(ndim,1),cof(lmax+1,nx,1)
C ... Local parameters
      integer i,ix,kkk

      do  10  i = 1, ndim
      ix = jp(i)
      do  10  kkk = 1, (lp1+1)*(lp2+1)
   10 cof(lxi+1,ix,kkk) = b(i,kkk)
      end
      subroutine on2phl(eph,lph,llink,rsm,el,add,r2,nr,idx,wk,phi,phil)
C- Make true radial parts of the linked phi's
      implicit none
C ... Passed parameters
      integer nr,lph,idx(nr,2),llink
      double precision wk(nr,0:1),phi(nr,0:1),phil(nr,0:1),
     .  r2(1),add(1),el,eph,rsm
C ... Local parameters
      integer ilm,l
      double precision tol
      parameter (tol=1d-10)

      call dpzero(phi,nr*(lph+1))
      call hansrg(rsm,lph,1,lph,eph,r2,nr,nr,tol,idx,wk,001,phi)
      if (llink == 0) return
      call hansrg(rsm,lph,1,lph,el,r2,nr,nr,tol,idx,wk,001,phil)
      do  10  l = 0, lph
      ilm = (l+1)**2
   10 call dpadd(phi(1,l),phil(1,l),1,nr,add(ilm))
      end
