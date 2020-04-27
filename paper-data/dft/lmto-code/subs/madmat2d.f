      subroutine madmat2D(nbas,npadl,npadr,s_lat,dmad)
C- Coefficients to Madelung matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,bas
Co Outputs
Co   dmad
Cr Remarks
Cr   Based on 2D-Electrostatics-Parry-SurfSci49,433.pdf
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,npadl,npadr
      double precision dmad(nbas,nbas)
C ... For structures
!       include 'structures.h'
      type(str_lat)::   s_lat
C ... Dynamically allocated local arrays
      real(8),allocatable :: pos(:,:)
C ... Local parameters
      integer ib,jb,i,stdo,ir1,ir,nbaspp
      double precision tau(3),dl,r1,ph,th,rotm(3,3),plat(3,3),qlat(3,3),vol,normal(3)
      real(8), parameter :: smallx=1d-8, pi=4*datan(1d0)
      procedure(integer) :: nglob,iprint
      procedure(real(8)) :: dlength

      call rx('not finished ... see 2D-Electrostatics-Parry-SurfSci49,433.pdf')

      nbaspp = nbas + 2*npadl + 2*npadr

C --- Rotate coordinate system so that plat(1)xplat(2) ~ qlat(3) on z axis ---
C      ph = datan2(s_lat%qlat(2,3),s_lat%qlat(1,3))
C      th = datan2(dlength(2,s_lat%qlat(1,3),1),s_lat%qlat(3,3))
C      call rotma(ph+pi/2,pi/2,th,rotm)
CC     call prmx('rotm',rotm,3,3,3)
C      call dgemm('N','N',3,3,3,1d0,rotm,3,s_lat%plat,3,0d0,plat,3)
CC     call prmx('rotated plat',plat,3,3,3)

C --- Unit vector normal to the plane
      call dinv33(s_lat%plat,1,qlat,vol)
      do  i = 1, 3
        normal(i) = qlat(i,3)/dlength(3,qlat(1,3),1)
      enddo


C ... Order positions by increasing z, z=distance to central plane

      allocate(pos(3,nbaspp))

C --- Generate Madelung matrix ---
      do  ib = 1, nbas
      do  jb = ib, nbas
        dl = 0
        forall (i=1:3) tau(i) = s_lat%pos(i,jb)-s_lat%pos(i,ib)
        call shortn(tau,tau,s_lat%dlv,s_lat%nkd)

C   ... Real space sums, backwards to sum small numbers
        if (tau(1)**2 + tau(2)**2 + tau(3)**2 > smallx) then
          ir1 = 1
        else
          ir1 = 2
          dl = dl - 2d0*s_lat%awald/dsqrt(pi)
        endif
        do  ir = s_lat%nkd, ir1, -1
          r1 = s_lat%alat*dsqrt((tau(1)-s_lat%dlv(1,ir))**2 +
     .                          (tau(2)-s_lat%dlv(2,ir))**2 +
     .                          (tau(3)-s_lat%dlv(3,ir))**2)
          dl = dl + derfc(s_lat%awald*r1)/r1
        enddo



        dmad(jb,ib) = dmad(ib,jb)

      enddo !jb
      enddo !ib

C --- Printout ---
      if (iprint() < 90) return
      stdo = nglob('stdo')
      write (stdo,1) nbas
    1 format(/' Madelung matrix: nbas=',i5)
      do  ib = 1, nbas
        print 2,(dmad(ib,jb),jb=1,nbas)
    2   format(7F11.7)
      enddo
      end
