      subroutine tetirr(qb,n1,n2,n3,ipq,ntet,idtet)
C-  Finds inequivalent tetrahedra and counts them
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ipq   :is a work array of dimension 6*n1*n2*n3.
Ci         :On input, ipq(i1,i2,i3) points to corresponding irreducible qp
Ci         :ipq is DESTROYED on output.
Ci   n1..n3:no. of divisions made along each reciprocal lattice vector
Ci   qb    :vectors of first microcell
Co Outputs:
Co   idtet :idtet(0,i) = number of tetrahedra of the i'th kind
Co         :idtet(1-4,i) points to the 4 irreducible k-points defining
Co         :the tetrahedron.
Co   ntet  :number of different tetrahedra
Cr Updates
Cu  02 Jun 13 completed replacement of f77 pointers
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer n1,n2,n3,ntet,ipq(n1,n2,n3),idtet(0:4,*)
      double precision qb(3,3)
C Local variables:
      integer i,j,i1,i2,i3,j1,j2,j3,k1,k2,k3,itet,mtet,ic,ii,ibtr(3,3),
     .  kcut(3,4,6),imc(0:1,0:1,0:1),iq(4),iprint
      double precision qb1(3,3)
      logical ipr
      integer, allocatable :: iwk(:)

c      call pshpr(40)

      ipr = iprint() >= 30 .and. 6*n1*n2*n3 > 1000
     .  .or. iprint() >= 40
      if (ipr) call info2(0,0,-1,' TETIRR: sorting %i tetrahedra ...',
     .  6*n1*n2*n3,0)

      call ccutup(qb,qb1,ibtr,kcut)
      ntet = 0
C --- Start loop over microcells ----
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
C   ... Set up identifiers at 8 corners of microcell
        do  k1 = 0, 1
          j1 = mod(i1+k1-1,n1) + 1
          do  k2 = 0, 1
            j2 = mod(i2+k2-1,n2) + 1
            do  k3 = 0, 1
              j3 = mod(i3+k3-1,n3) + 1
              imc(k1,k2,k3) = ipq(j1,j2,j3)
            enddo
          enddo
        enddo

C   --- Start loop over tetrahedra ---
        do  itet = 1, 6
          do  ic = 1, 4
            k1 = kcut(1,ic,itet)
            k2 = kcut(2,ic,itet)
            k3 = kcut(3,ic,itet)
            iq(ic) = imc(k1,k2,k3)
          enddo

C    ... Order the identifiers
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
          idtet(0,ntet) = 1
          do  i = 1, 4
            idtet(i,ntet) = iq(i)
          enddo
        enddo
      enddo
      enddo
      enddo

C --- Eliminate duplicate tetrahedra ---
      mtet = ntet
C ... Sort it with fast ivheap unless there is not enough memory
      if (.true.) then
        call ivheap(5,mtet,idtet,ipq,1)
        allocate(iwk(mtet*5))
        call ivprm(5,mtet,idtet,iwk,ipq,1)
        deallocate(iwk)
*       call yprm('idtet',0,idtet,0,5,5,mtet)
      else
        call ivshel(5,mtet,idtet,ipq,.false.)
*       call yprm('idtet',0,idtet,0,5,5,mtet)
      endif
      ntet = 1
      do  i = 2, mtet
C        print 337, i, idtet(1,i), idtet(2,i), idtet(3,i), idtet(4,i)
C  337   format(i4,4i6)
        if (idtet(1,i) == idtet(1,ntet)  .and.
     .      idtet(2,i) == idtet(2,ntet)  .and.
     .      idtet(3,i) == idtet(3,ntet)  .and.
     .      idtet(4,i) == idtet(4,ntet)) then
          idtet(0,ntet) = idtet(0,ntet)+1
        else
          ntet = ntet+1
          idtet(1,ntet) = idtet(1,i)
          idtet(2,ntet) = idtet(2,i)
          idtet(3,ntet) = idtet(3,i)
          idtet(4,ntet) = idtet(4,i)
        endif
      enddo

      if (ipr) call info2(0,0,0,' %i inequivalent ones found',ntet,0)
      end
