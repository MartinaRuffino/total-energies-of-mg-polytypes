      subroutine tetfbz(qlat,n1,n2,n3,nadd,skipgammacell,ngcell,half,defrm,
     .  idtet,qbzw,ib1bz,ntet)
C- Finds tetrahedra in entire 1st BZ, without symmetry operations
C ----------------------------------------------------------------------
Ci Inputs:
Ci  qlat     :reciprocal lattice vectors
Ci  n1,n2,n3 :number of divisions along the three R.L.V.
Ci  lshft    :For each recip. lattice vector :
Ci           :0 center the mesh through the origin (k=0)
Ci           :1 center the mesh straddling the origin --- no point at 0.
Ci  nq, no. of irreducible k-points;
Co Outputs:
Co   idtet :idtet(1:4,i) points to the 4 k-points defining the ith
Co         :tetrahedron.  Note: idtet refers to the padded qbzw
Co   qbzw  :points on the regular mesh, padded to repeat points
Co         :on the boundary lines
Co   ib1bz :maps the points in the padded qp list to the original list
Cr Remarks
Cr   tetfbz is designed as setup for routines that calculate weights
Cr   for response functions, e.g.
Cr     \int d^3k f(e(k)) (1-f(e(q+k))) \delta (omega - e(q+k) + e(k) )
Cr   where f(E) is the Fermi function.
Cr   It is related to tetirr.f, but has the following differences:
Cr     1. No symmetry operations are assumed (q may be aribrary)
Cr     2. k-points are returned in qbzw.  qbzw is enalarged to include
Cr     boundary points twice (faciliting treatment of periodicity)
Cr     3. index ib1bz is returned.
Cr     4. idtet is dimensioned differently
Cu Updates
Cu   27 Apr 19 Extra arguments in preparation for more new GW
Cu   04 Sep 09 Adapted from Takao's tetfbzf
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical skipgammacell
      integer n1,n2,n3,nadd,ngcell
      integer ib1bz((n1+1)*(n2+1)*(n3+1)),idtet(4,6*n1*n2*n3)
      double precision defrm,qlat(3,3),qbzw(3,(n1+1)*(n2+1)*(n3+1)),half(3)
C ... Local parameters
      integer i,i1,i2,i3,ic,ii,itet,j,k1,k2,k3,kount,ntet,
     .  imc(0:1,0:1,0:1),kcut(3,4,6),ibtr(3,3),iq(4)
      integer indexkw(0:n1,0:n2,0:n3),indexk(0:n1-1,0:n2-1,0:n3-1)
      double precision QB(3,3),QB1(3,3)
C     double precision qbzx(3)
C     double precision epss
C     parameter (epss = 1d-12)
C     logical,save ::chk=.true.

      print *, 'tetfbz',half

      if (nadd /= 0 .and. nadd /= 1) call rx('tetfbz : illegal nadd')
      if (defrm /= 1) call rx('tetfbz : not ready for defrm')
      if (skipgammacell) call rx('tetfbz : not ready for skipgammacell')

C ... Offsets for each of 3 directions
C      hf(1) = 0d0
C      hf(2) = 0d0
C      hf(3) = 0d0
C      if (lshft(1) /= 0) hf(1) = 0.5d0
C      if (lshft(2) /= 0) hf(2) = 0.5d0
C      if (lshft(3) /= 0) hf(3) = 0.5d0

C ... Offsets for each of 3 directions
      qb(1:3,1) = qlat(1:3,1)/n1
      qb(1:3,2) = qlat(1:3,2)/n2
      qb(1:3,3) = qlat(1:3,3)/n3

C ... Index to long vector of kpts in 1st BZ.
      kount = 0
      do  i1 = 1, n1+nadd
        do  i2 = 1, n2+nadd
          do  i3 = 1, n3+nadd
            kount = kount + 1
            indexk(i1-1,i2-1,i3-1) = kount
C            qbzx(1:3)= qb(1:3,1)*(i1-1+half(1)) +
C     .                 qb(1:3,2)*(i2-1+half(2)) +
C     .                 qb(1:3,3)*(i3-1+half(3))
C            write(*,"(3i4,3d26.18)") i1,i2,i3,qbzx
          end do
        end do
      end do

C ... Padded vector of kpts in 1st BZ, including both boundary points
      kount = 0
      do  i1 = 1, n1+1
        do  i2 = 1, n2+1
          do  i3 = 1, n3+1
            kount = kount + 1
            indexkw(i1-1,i2-1,i3-1) = kount
            qbzw(1:3,kount) =
     .        qb(1:3,1)*(i1-1+half(1)) +
     .        qb(1:3,2)*(i2-1+half(2)) +
     .        qb(1:3,3)*(i3-1+half(3))
            if(nadd==0) ib1bz(kount) = indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
            if(nadd==1) ib1bz(kount) = indexk(i1-1, i2-1, i3-1)
          end do
        end do
      end do

      call ccutup(qb,qb1,ibtr,kcut)
      ntet = 0
C --- Start loop over microcells ----
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
C   ... Set up identifiers at 8 corners of microcell
C        do  k1 = 0, 1
C          j1 = i1 -1 + k1
C          do  k2 = 0, 1
C            j2 = i2 -1 + k2
C            do  k3 = 0, 1
C              j3 = i3 -1 + k3
C              imc(k1,k2,k3) = indexkw(j1,j2,j3)
C            enddo
C          enddo
C        enddo
        imc(0:1,0:1,0:1) = indexkw(I1-1:I1,I2-1:I2,I3-1:I3)

C --- Start loop over tetrahedra ---
        do  itet = 1, 6
          do  ic = 1, 4
            k1 = kcut(1,ic,itet)
            k2 = kcut(2,ic,itet)
            k3 = kcut(3,ic,itet)
            iq(ic) = imc(k1,k2,k3)
          enddo

C  ... Order the identifiers
          do  j = 1, 3
            do  i = 1, 4-j
              if (iq(i) > iq(i+1)) then
                ii = iq(i)
                iq(i) = iq(i+1)
                iq(i+1) = ii
              endif
            enddo
          enddo
          ntet = ntet+1
C         idtet(0,ntet) = 1
          do  i = 1, 4
            idtet(i,ntet) = iq(i)
          enddo
        enddo
      enddo
      enddo
      enddo

C      do  i = 1, ntet
C        print 337, i, idtet(1,i), idtet(2,i), idtet(3,i), idtet(4,i)
C  337   format(i4,4i6)
C      enddo
C      print "(1x,' tetfbz: ',2i8 )",  ntet, 6*n1*n2*n3
      END

C      subroutine tetfbz(qlat,n1,n2,n3,lshft,idtet,qbzw,ib1bz)
CC- Finds tetrahedra in entire 1st BZ, without symmetry operations
CC ----------------------------------------------------------------------
CCi Inputs:
CCi  qlat     :reciprocal lattice vectors
CCi  n1,n2,n3 :number of divisions along the three R.L.V.
CCi  lshft    :For each recip. lattice vector :
CCi           :0 center the mesh through the origin (k=0)
CCi           :1 center the mesh straddling the origin --- no point at 0.
CCi  nq, no. of irreducible k-points;
CCo Outputs:
CCo   idtet :idtet(1:4,i) points to the 4 k-points defining the ith
CCo         :tetrahedron.  Note: idtet refers to the padded qbzw
CCo   qbzw  :points on the regular mesh, padded to repeat points
CCo         :on the boundary lines
CCo   ib1bz :maps the points in the padded qp list to the original list
CCr Remarks
CCr   tetfbz is designed as setup for routines that calculate weights
CCr   for response functions, e.g.
CCr     \int d^3k f(e(k)) (1-f(e(q+k))) \delta (omega - e(q+k) + e(k) )
CCr   where f(E) is the Fermi function.
CCr   It is related to tetirr.f, but has the following differences:
CCr     1. No symmetry operations are assumed (q may be aribrary)
CCr     2. k-points are returned in qbzw.  qbzw is enalarged to include
CCr     boundary points twice (faciliting treatment of periodicity)
CCr     3. index ib1bz is returned.
CCr     4. idtet is dimensioned differently
CCu Updates
CCu   04 Sep 09 Adapted from Takao's tetfbzf
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer n1,n2,n3,lshft(3)
C      integer ib1bz((n1+1)*(n2+1)*(n3+1)),idtet(4,6*n1*n2*n3)
C      double precision qlat(3,3),qbzw(3,(n1+1)*(n2+1)*(n3+1))
CC ... Local parameters
C      integer i,i1,i2,i3,ic,ii,itet,j,j1,j2,j3,k1,k2,k3,kount,ntet,
C     .  imc(0:1,0:1,0:1),kcut(3,4,6),ibtr(3,3),iq(4)
C      integer indexkw(0:n1,0:n2,0:n3),indexk(0:n1-1,0:n2-1,0:n3-1)
C      double precision QB(3,3),QB1(3,3),hf(3)
CC     double precision qbzx(3)
CC     double precision epss
CC     parameter (epss = 1d-12)
CC     logical,save ::chk=.true.
C
CC ... Offsets for each of 3 directions
C      hf(1) = 0d0
C      hf(2) = 0d0
C      hf(3) = 0d0
C      if (lshft(1) /= 0) hf(1) = 0.5d0
C      if (lshft(2) /= 0) hf(2) = 0.5d0
C      if (lshft(3) /= 0) hf(3) = 0.5d0
C
CC ... Offsets for each of 3 directions
C      qb(1:3,1) = qlat(1:3,1)/n1
C      qb(1:3,2) = qlat(1:3,2)/n2
C      qb(1:3,3) = qlat(1:3,3)/n3
C
CC ... Index to long vector of kpts in 1st BZ.
C      kount = 0
C      do  i1 = 1, n1
C        do  i2 = 1, n2
C          do  i3 = 1, n3
C            kount = kount + 1
C            indexk(i1-1,i2-1,i3-1) = kount
CC            qbzx(1:3)= qb(1:3,1)*(i1-1+hf(1)) +
CC     .                 qb(1:3,2)*(i2-1+hf(2)) +
CC     .                 qb(1:3,3)*(i3-1+hf(3))
CC            write(*,"(3i4,3d26.18)") i1,i2,i3,qbzx
C          end do
C        end do
C      end do
C
Ccccccccccccccccccccccccccccc
Cc      DO  j1 = 0,n1-1
Cc      DO  j2 = 0,n2-1
Cc      DO  j3 = 0,n3-1
Cc        write(6,"(' j1j2j3=',3i4,' ix=',i6)") J1,J2,J3,indexk(J1,J2,J3)
Cc      enddo
Cc      enddo
Cc      enddo
Cc      stop
Ccccccccccccccccccccccccccccc
C
CC ... Padded vector of kpts in 1st BZ, including both boundary points
C      kount = 0
C      do  i1 = 1, n1+1
C        do  i2 = 1, n2+1
C          do  i3 = 1, n3+1
C            kount  = kount + 1
C            indexkw(i1-1,i2-1,i3-1) = kount
C            qbzw(1:3,kount) =
C     .        qb(1:3,1)*(i1-1+hf(1)) +
C     .        qb(1:3,2)*(i2-1+hf(2)) +
C     .        qb(1:3,3)*(i3-1+hf(3))
C            ib1bz(kount) =
C     .        indexk(mod(i1-1,n1), mod(i2-1,n2), mod(i3-1,n3))
C          end do
C        end do
C      end do
C
C      call ccutup(qb,qb1,ibtr,kcut)
C      ntet = 0
CC --- Start loop over microcells ----
C      do  i3 = 1, n3
C      do  i2 = 1, n2
C      do  i1 = 1, n1
CC   ... Set up identifiers at 8 corners of microcell
C        do  k1 = 0, 1
C          j1 = i1 -1 + k1
C          do  k2 = 0, 1
C            j2 = i2 -1 + k2
C            do  k3 = 0, 1
C              j3 = i3 -1 + k3
C              imc(k1,k2,k3) = indexkw(j1,j2,j3)
C            enddo
C          enddo
C        enddo
C
CC --- Start loop over tetrahedra ---
C        do  itet = 1, 6
C          do  ic = 1, 4
C            k1 = kcut(1,ic,itet)
C            k2 = kcut(2,ic,itet)
C            k3 = kcut(3,ic,itet)
C            iq(ic) = imc(k1,k2,k3)
C          enddo
C
CC  ... Order the identifiers
C          do  j = 1, 3
C            do  i = 1, 4-j
C              if (iq(i) > iq(i+1)) then
C                ii = iq(i)
C                iq(i) = iq(i+1)
C                iq(i+1) = ii
C              endif
C            enddo
C          enddo
C          ntet = ntet+1
CC         idtet(0,ntet) = 1
C          do  i = 1, 4
C            idtet(i,ntet) = iq(i)
C          enddo
C        enddo
C      enddo
C      enddo
C      enddo
C
CC      do  i = 1, ntet
CC        print 337, i, idtet(1,i), idtet(2,i), idtet(3,i), idtet(4,i)
CC  337   format(i4,4i6)
CC      enddo
CC      print "(1x,' tetfbz: ',2i8 )",  ntet, 6*n1*n2*n3
C      END
