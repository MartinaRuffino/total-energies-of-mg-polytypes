      subroutine stonerpb(nq,nw,nmbas,nbloch,
     .  qp,momsite,mmnorm,emesh,zxq)
C=============================================================
C ---  <~e| X | e~>  matrix calculation. by using X^{-1}=X0^{-1}+U
C ----------------------------------------------------------------------
Ci Inputs
Ci   nq       : total number of q points in calculation of X(q,w)
Ci   nw       : total number of w points in calculation of X(q,w)
Ci   nmbas    : number of magnetic sites  = nmag
Ci   qp       : qp(3,nq)  vector of each q (not used here)
Ci   momsite  : magnetic moment m (not used here)
Ci   mmnorm   : <m|m> (not used here)
Ci   eiqrm    : eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
Ci   emesh    : energy w mesh
Ci   x0et      : <~e|X0|~e>
Co Outputs:
Cl Local variables
Cl   rho0    : radial part of true electron density on each mesh (xyz) = rho1/r**2
Cl   rho0xyz : true electron density on each mesh point specified by (ir,ip)
Cl             rho0xyz(ir,ip,isp)=Sum_{ilm}{ rho0(ir,ilm,isp)*yl(ip,ilm) }
Cl   vxc     : vxc(ir,ip,1:2) for spin_1 and spin_2 ;
Cl             vxc(:,:,4) = ( vxc(:,:,1)-vxc(:,:,2) ) / 2
Cl             vxc(:,:,3) = idol
Cl   bom     : bom(ir,ip) = b(r)/m(r)
Cl   sum_bibj: <B_i | B_j> integrate inside sphere
Cl   sum_bibombj: <B_i | b(r)/m(r) | B_j>
Cr Remarks
Cr    In ASA, m(r_) is a radial function, but not in FP.
Cr    we need to build a spherical mesh, and sum the rho on each mesh
Cr    point.
Cr    Vxc generated from LMTO package is in units of Ryberg.
Cr    Numerical spherical mesh is tabulated in wxp.chk file, which is generated in LMTO
Cr    package ./lmfsph
Cu Updates
Cu    20 May 08 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nq,nw,nmbas,nbloch
      real(8) momsite(nmbas), mmnorm(nmbas),emesh(nw),qp(3,nq)
      complex(8) zxq(nbloch,nbloch,nq,nw)
C ... Local parameters
      integer jb,iw,iq
      real(8) freq(nw)

C     nwx = 400

c      allocate( x0inv(nmbas,nmbas,nq,nwx) )
c      allocate( xinv(nmbas,nmbas,nq,nwx), x2et(nmbas,nmbas,nq,nwx) )
c      allocate( xinvh(nmbas,nmbas,nq,nwx),mxevl_xinvh(nq,nwx) )
c      allocate( dxidw(nmbas,nmbas,nq,nwx))
c      allocate( dxidw_eb(nmbas,nmbas,nq,nwx))
c      allocate( mxevl2_xinvh(nq,nwx) )

C ... Sanity check
      if ( emesh(1) /= 0d0 ) call rx('stonerpb: w(1) /= 0')
      freq=emesh

c./lmf lsmo56 '--chimedit~new 2~read tkrs'

      jb=1
      open(111, file='X0pbqw.allqb')
      write(111,"( 4i4 )") nq,nmbas,nw
      do  iw = 1, nw
        write(111,301)  (zxq( jb , jb,iq,iw), iq=1,nq )
      enddo
  301 format(1000d23.15)
      end subroutine stonerpb
