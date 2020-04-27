      subroutine plch1(ip,vx,vy,x1,x2,y1,y2,nx,ny,pos0,dist,np,x,y,z,a)
C- set up list of points in a plane for chden contour plot.
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs: vx, vy, vz cartesian set in which vx is the direction
Co          between atoms 1 and 2 in the list of three atoms that
Co          define the plane of plotting. vy lies in the xy plane
Co          defined by the three atoms.
Co          A(i,j) is the transformation matrix for writing the
Co          atomic coordinates in terms of the cartesian set vx, vy, vz
Co
Cr Remarks
C ----------------------------------------------------------------------
      implicit real*8 (a-h,p-z), integer(o)
      dimension vx(3),vy(3),vz(3),x(nx,ny),y(nx,ny),z(nx,ny),pos0(3),
     .   xx(3),a(3,3)
c ---- orthonormalize vectors which determine the plane ----
      call dpdot(vx,vy,3,axy)
      call dpdot(vx,vx,3,axx)
      do 1 m=1,3
  1   vy(m)=vy(m)-vx(m)*(axy/axx)
      call dpdot(vx,vx,3,axx)
      call dpdot(vy,vy,3,ayy)
      do 2 m=1,3
      vx(m)=vx(m)/dsqrt(axx)
  2   vy(m)=vy(m)/dsqrt(ayy)
      call cross(vx,vy,vz)
c ---- permute according to ip ----------
      i1=ip/100
      mm=ip-i1*100
      i2=mm/10
      i3=mm-i2*10
      write(6,723) i1,i2,i3
  723 format(/' plch1:  vec permutation indices:',3i4)
      if(max0(i1,i2,i3) > 3.or.min0(i1,i2,i3) < 1)
     .    call rx('plch1: vec permutation index not 1,2 or 3')
      do 40 m=1,3
      xx(1)=vx(m)
      xx(2)=vy(m)
      xx(3)=vz(m)
      vx(m)=xx(i1)
      vy(m)=xx(i2)
  40  vz(m)=xx(i3)
      write(6,601) 'vx',vx
      write(6,601) 'vy',vy
      write(6,601) 'vz',vz
  601 format(' vector ',a2,'=',3f12.6)
c ------ make arrays x,y,z ---------
      np=nx*ny
      do 11 i=1,nx
      x0=(i-1)*(x2-x1)/(nx-1d0)+x1
      do 11 j=1,ny
      y0=(j-1)*(y2-y1)/(ny-1d0)+y1
      x(i,j)=x0*vx(1)+y0*vy(1)+dist*vz(1)+pos0(1)
      y(i,j)=x0*vx(2)+y0*vy(2)+dist*vz(2)+pos0(2)
  11  z(i,j)=x0*vx(3)+y0*vy(3)+dist*vz(3)+pos0(3)
C --- Make tranformation matrix ---
      a(1,1) = vx(1)
      a(1,2) = vy(1)
      a(1,3) = vz(1)
      a(2,1) = vx(2)
      a(2,2) = vy(2)
      a(2,3) = vz(2)
      a(3,1) = vx(3)
      a(3,2) = vy(3)
      a(3,3) = vz(3)
      write(*,602)((a(i,j),i=1,3),j=1,3)
  602 format (' Transformation matrix: ',3f10.6,2(/24x,3f10.6))
      end
