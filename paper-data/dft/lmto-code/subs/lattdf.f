      subroutine lattdf(ldist,defgrd,plat,nbas,bas,ngen,gen,ag)
C- Rotates or deforms lattice
C ----------------------------------------------------------------
Ci Inputs
Ci   ldist: 0, no deformations.  For abs(ldist):
Ci          1: defgrd holds rot about spec'd angle
Ci          2, lattice deformed with a general linear transformation
Ci          3, lattice deformed by a shear.
Ci          SIGN ldist <0 => suppress rotation of plat
Ci   nbas  :size of basis
Ci   ngen  :number of generators of the symmetry group
Co Outputs
Cio  defgrd:transformation matrix, whose form is specified by ldist
Cio        :If ldist ==1 or ldist=3, defgrd is returned as a matrix
Cio  plat  :primitive lattice vectors, in units of alat
Cio        :On output, plat is transformed by defgrd
Cio  bas   :basis vectors, in units of alat
Cio        :On output, bas is transformed by defgrd
Cio  gen   :On input, generators of the symmetry group
Cio        :On output, generators are transformed by defgrd
Cio        :Note: gen can only be transformed by a rotation
Cio  defgrd:Output defgrd is actual transformation matrix.
Cu Updates
Cu   01 Jun 15 Also rotate translation part of symmetry group generators
Cu   19 Mar 06 Blend Voigt strains into other types
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldist,nbas,ngen
      double precision defgrd(3,3),plat(3,3),bas(3,1),gen(3,3,1),ag(3,1)
C ... Local parameters
      double precision work(3,3),rinv(3,3),det,gold(3,3),agold(3)
      integer ipr,j,ib,stdo,nglob
C ... External calls
      external dcopy,defrr,dinv33,dmpy,getpr,makrot,shear
C     Debugging
C      integer nlx,nl2
C      parameter (nlx=4,nl2=nlx*nlx)
C      double precision rYL(nl2,nl2)

      if (ldist == 0) return

      stdo = nglob('stdo')
      call getpr(ipr)
      call dcopy(9,defgrd,1,work,1)
      if (iabs(ldist) == 1) call makrot(work,defgrd)
      if (iabs(ldist) == 3) then
        det = defgrd(1,3)
        if (det == 0) return
        call shear(0,plat,bas,det,defgrd,defgrd)
      endif
      call info2(30,1,0,' LATTDF:  deformation matrix for mode %i',ldist,0)
      do  j = 1, 3
        call info2(30,0,0,'%3;12,7D',defgrd(1,j),0)
      enddo

C     call prmx('rotation matrix',defgrd,3,3,3)
C     call ylmrtg(nl2,defgrd,rYL)
C     call prmx('rothrm: rYL',rYL,nl2,nl2,nl2)

C ... Rotate or shear plat
      if (ldist >= 1 .and. ldist <= 3) then
        call dcopy(9,plat,1,work,1)
        call dinv33(defgrd,0,rinv,det)
        call dmpy(defgrd,3,1,work,3,1,plat,3,1,3,3,3)
        call info0(30,0,0,'%11pLattice vectors:%25fTransformed to:')
        do  j = 1, 3
          call info2(30,0,0,'%3;12,7D  %3;12,7D',work(1,j),plat(1,j))
        enddo
      endif

      if (iabs(ldist) >= 1 .and. iabs(ldist) <= 3) then
        if (nbas > 0)
     .  call info0(30,0,0,'%12pBasis vectors:%26fTransformed to:')
        do  ib = 1, nbas
          call dcopy(3,bas(1,ib),1,work,1)
          call dmpy(defgrd,3,1,work,3,1,bas(1,ib),3,1,3,1,3)
C          if (ipr >= 30) print '(3f12.7,2x,3f12.7)',
C     .      (work(i,1), i=1,3), (bas(i,ib), i=1,3)
          call info2(30,0,0,'%3;12,7D  %3;12,7D',work(1,1),bas(1,ib))
        enddo

        if (ngen > 0) call info0(30,1,0,'%15pGroup op:%28fRotated to:')
        call dinv33(defgrd,0,rinv,det)
        do  ib = 1, ngen
          call dcopy(9,gen(1,1,ib),1,gold,1)
          call dmpy(defgrd,3,1,gen(1,1,ib),3,1,work,3,1,3,3,3)
          call dmpy(work,3,1,rinv,3,1,gen(1,1,ib),3,1,3,3,3)
C          if (ipr >= 30) print '(/(3f12.7,2x,3f12.7))',
C     .      ((gold(j,i), i=1,3), (gen(j,i,ib), i=1,3),j=1,3)
          do  j = 1, 3
            call info2(30,0,0,'%3;12,7D  %3;12,7D',gold(1,j),gen(j,:,ib))
          enddo
          call dcopy(3,ag(1,ib),1,agold,1)
          call dmpy(defgrd,3,1,agold,3,1,ag(1,ib),3,1,3,1,3)
          call info2(30,0,1,'%3;12,7D  %3;12,7D   ag',agold,ag(1,ib))
        enddo

      endif
      end
C      subroutine rotp(rotm,np,p,prot)
CC- Rotates a vector of points by a rotation matrix
CC  Returns prot = rotm*p.  prot may point to the same address space as p
C      implicit none
C      integer np
C      double precision p(3,np),prot(3,np),rotm(3,3),h(3)
C      integer i,j,ip
C
C      do  10  ip = 1, np
C        do  1  i = 1, 3
C        h(i) = 0
C        do  1  j = 1, 3
C    1   h(i) = h(i) + rotm(i,j)*p(j,ip)
C        do  2  i = 1, 3
C    2   prot(i,ip) = h(i)
C   10 continue
C      end
