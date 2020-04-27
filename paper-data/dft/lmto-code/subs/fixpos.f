      subroutine fixpos(pos,nbas,tol,ng,plat,g,ag,istab)
C- Adjusts site positions to agree with given symmmetry to machine precision
C ----------------------------------------------------------------------
Ci Inputs:
Ci   pos:   basis vectors (scaled by alat)
Ci   nbas:  number of atoms in the basis
Ci   tol:   largest separation at which atoms are considered coincident
Ci   ng:    number of symmetry operations
Ci   plat:  primitive lattice vectors (scaled by alat)
Ci   g,ag:  point and translation group operators
Ci   istab: atom transformation table; see symtab
Co Outputs:
Co   pos:   basis vectors are adjusted.
Cr Remarks:
Cr   Generally atoms of the same class do not sit exactly on
Cr   symmetry-related positions. In this subroutine each atomic
Cr   position is replaced by the average of the position itself and
Cr   the generated positions of the atoms of the same class.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters:
      integer nbas,ng,istab(nbas,ng)
      double precision pos(3,*),plat(*),g(9,*),ag(3,*),tol
C Local parameters:
      integer ibas,jbas,m,ig,lgunit
      double precision dbas(3),bast(3),sum,tol2,qlat(3,3),vol,ddot
      double precision sdpos(3,nbas)

      tol2 = 2*tol
      call dpzero(sdpos,3*nbas)
      sum = 0
      call dinv33(plat,1,qlat,vol)
      do  ibas = 1, nbas
        do  ig = 1, ng
          jbas = istab(ibas,ig)
          call dmpy(g(1,ig),3,1,pos(1,ibas),3,1,bast,3,1,3,1,3)
          do  m = 1, 3
            dbas(m) = bast(m)+ag(m,ig) - pos(m,jbas)
          enddo
C         print 333, ' input basj',pos(1,jbas),pos(2,jbas),pos(3,jbas)
C         print 333, ' input dbas',dbas
    1     format(a,3F12.6)
          call shorbz(dbas,dbas,plat,qlat)
c         print 334, 'output dbas', dbas
C     ... Debugging check
          sum = sum + abs(dbas(1))+abs(dbas(2))+abs(dbas(3))
          if (abs(dbas(1)) > tol2 .or. abs(dbas(2)) > tol2.or.
     .        abs(dbas(3)) > tol2) call fexit(-1,111,
     .      'Exit -1 FIXPOS: positions incompatible with symgrp:'//
     .      '  dpos=%d',max(dbas(1),dbas(2),dbas(3)))
          if (abs(dbas(1)) > tol .or. abs(dbas(2)) > tol .or.
     .        abs(dbas(3)) > tol) call awrit4(
     .      ' FIXPOS (warning): sites %i,%i incompatible '//
     .        'with grp op %i:  dpos=%d',' ',80,lgunit(1),ibas,jbas,ig,
     .        max(dbas(1),dbas(2),dbas(3)))
    2     format(a,3F18.12)
          call daxpy(3,1d0,dbas,1,sdpos(1,jbas),1)
        enddo
      enddo
      sum = dsqrt(ddot(3*nbas,sdpos,1,sdpos,1)/3/nbas)
      call daxpy(3*nbas,1d0/ng,sdpos,1,pos,1)
      call info2(20,0,0,' FIXPOS: shifted site positions by average %;3g',sum/ng,0)

C     call prmx('pos',pos,3,3,nbas)
C     pause

      end
