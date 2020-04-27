      subroutine stewld(as,tol,alat0,alat,plat,gt,lmax,ixp,dxp)
C- Set up for Ewald summations (molecules version)
C  Scalar data passed through arrays ixp and dxp.  Conventions:
C      dxp       ixp
C 1:   alat      flag set to 1 indicating nonzero on-site terms
      implicit none
      integer ixp(5),lmax
      double precision alat0,alat,as,plat(3,3),dxp(5),gt,tol
      double precision rb(3,3),qb(3,3),qlat(3,3),qdist0,rdist0,
     .  tripl,xx,a0,tol1,awald,g1,g2,g3,r0,radd,q0,qadd,vol,vol0
      parameter (g1=1d0,g2=1d0,g3=1d0)
      integer odlat,orlat,owk,ipr,k,m,nkd,nkr,nkdest,nkrest
      real w(1)
      common /w/ w

      alat = 1d0
      dxp(1) = alat
      ixp(1) = 0

      end
