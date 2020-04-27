      subroutine wsloc(wscrk,nn,v1,v2,volbz,qp,nkp,nkxyz,wtkp,wscrl)
C- Make local part of ASA SX potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   wscrk :k-dependent screened exchange
Ci   nn    :dimension of W (nbas or sum 1+lmax)
Ci   v1    :
Ci   v2    :
Ci   volbz :BZ volume
Ci   qp    :k-point
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nkxyz :no of divisions made along each reciprocal lattice vector
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Co Outputs
Co   wscrl :local part of W
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer nn,nkp,nkxyz(3)
      double precision wscrk(nn,nn,2,nkp),wscrl(nn,nn,2),
     .  v1(nn,nn,2),v2(nn,nn,2),qp(3,*),wtkp(*),volbz
      double precision fpi,q0,q2,w,dsum,v,wg
      integer i,j,iq,lgunit

      fpi = 16*datan(1d0)
      q0 = (volbz*3d0/fpi/nkxyz(3)/nkxyz(2)/nkxyz(1))**(1d0/3)
      q2 = dsqrt(qp(1,2)**2 + qp(2,2)**2 + qp(3,2)**2)
      print '(/a,2f10.4)', ' wsloc: q0, q2:        ',q0,q2

C     call yprm('wstat(q=0)',2,wscrk,nn*nn,nn,nn,nn)
C     call yprm('wstat(q=q2)',2,wscrk(1,1,1,2),nn*nn,nn,nn,nn)
C     call yprm('v1',2,v1,nn*nn,nn,nn,nn)

C --- Est. epsilon as sum_nn (w(iq=2)-w(iq=1))/(v1(iq=2)-v1(iq=1)) ---
      w = (dsum(nn,wscrk(1,1,1,2),nn+1) - dsum(nn,wscrk,nn+1))/nn
      v = (dsum(nn,v2,nn+1) - dsum(nn,v1,nn+1))/nn
      print '(a,3f10.4)', ' wsloc:  w, v, epsilon:', w,v,v/w
      i = lgunit(2)
      write (i,'('' eps'', f8.2)') v/w

C --- Replace w(q=0) ---
C     Is this to estimate the part left out because the 1/q**2
C     contribution to W was suppressed??
      wg = 3*w*(q2/q0)**2
      wg = wg/2
      do  20  i = 1, nn
      do  20  j = 1, nn
   20 wscrk(j,i,1,1) = wscrk(j,i,1,1) + wg

C --- Accumulate W_RR' ;  See Rucker's preprint Eq. 16. ---
      call dpzero(wscrl,nn**2*2)
      do  30  iq = 1, nkp
        w = abs(wtkp(iq))/2d0
        do  32  i = 1, nn
        do  32  j = 1, nn
        wscrl(j,i,1) = wscrl(j,i,1) + wscrk(j,i,1,iq)*w
   32   wscrl(j,i,2) = wscrl(j,i,2) + wscrk(j,i,2,iq)*w
   30 continue

      end








