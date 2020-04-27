      subroutine vesft(s_lat,ng,gv,kv,k1,k2,k3,smrho,cv,sum)
C- Make electrostatic potential of density given in recip space
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
Ci   smrho :FT of smooth density on uniform mesh
Co Outputs
Co   cv    :FT of smooth electrostatic potential in glist form,
Co         :is added to cv
Co   sum   :integral density*(estat pot + input cv)
Cr Remarks
Cu Updates
Cu   16 Jun 16 estat potential is now added to smpot, and stored in cv
Cu             To reproduce old results, initialize input cv to zero.
Cu   10 Nov 11 Begin migration to f90 structures
Cu   22 Apr 00 Adapted from nfp ves_ft.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ng
      integer kv(ng,3),k1,k2,k3
      double precision gv(ng,3),sum
      double complex smrho(k1,k2,k3),cv(ng)
C ... For structures
!      include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      complex(8), allocatable :: cn(:)
C ... Local parameters
      integer i
      double precision pi,pi8,alat,vol,tpiba,g2
      double complex ccc

      call tcn('vesft')
      alat = s_lat%alat
      vol = s_lat%vol
      pi   = 4d0*datan(1d0)
      pi8  = 8*pi
      tpiba=2*pi/alat

C ... Gather density coefficients
      allocate(cn(ng))
      call gvgetf(ng,1,kv,k1,k2,k3,smrho,cn)
C     call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv0)

C ... smpot(G) = 8 pi /G**2 smrho(G)
C     call info2(30,1,0,' vesft:  smooth density coeff to (0,0,0) '//
C    .  'is %;12,6D',cv(1),0)
      sum = 0d0
      do  i = 2, ng
        g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
        ccc = (pi8/g2)*cn(i) + cv(i)
        sum = sum + dconjg(cn(i))*ccc
        cv(i) = ccc
      enddo
      sum = vol*sum

      deallocate(cn)
      call tcx('vesft')
      end
