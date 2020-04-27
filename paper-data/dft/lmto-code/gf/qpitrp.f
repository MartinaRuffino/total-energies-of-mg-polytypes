      subroutine qpitrp(n1,n2,n3,ifac,qb,qp,f,fqp)
C- Linearly interpolate a function on tabulated on a mesh.
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1..n3:number of mesh points along qlat vectors
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          (not used here)
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   qp    :k-point at which to evaluate f
Co Outputs
Co   fqp   :result of f at qp
Cr Remarks
Cr   Algorithm:
Cr  *qp is mapped into a multiple of qb, the primitive mesh vectors
Cr  *Mesh point closest to point qp is found; call this point q4.
Cr  *Find q1,q2,q3 = three nearest neighbors of the cube, i.e.
Cr   q4 +/- (1,0,0), (0,1,0), (0,0,1)
Cr  *Find f1,f2,f3 = function values at points q1,q2,q3
Cr  *Make array q_mi = ith component of vector qm - q4
Cr  *Calculate function value at qp as
Cr     f(qp) = (qp-q4)_i (q_mi)^-1 f_m + f4
Cr   NB: qmi is diagonal, with elements +/- 1, so
Cr     f(qp) = (qp-q4)_i (q_ii) f_i + f4
Cr  Debugging:
Cr    mc plat -i -t >qlat
Cr    mc -vn=12 qlat -s1/n >qb
Cu Updates
Cu   24 Feb 07 Bug fix for noncubic lattices
C ----------------------------------------------------------------------
      implicit none
      integer n1,n2,n3,ifac(3)
      double precision qp(3),qb(3,3),f(n1,n2,n3),fqp
      integer nk(3),i,ii(4,3),m
      double precision qbi(3,3),det,qqb(3),xl(3),fm(4),qii(3),qpmq4(3)
C     integer jj1,jj2,jj3,k
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
C     .                    (jj2*ifac(2)-1)*qb(k,2) +
C     .                    (jj3*ifac(3)-1)*qb(k,3)


C      qp(1) = .62d0
C      qp(2) = .22d0
C      qp(3) = .14d0
C      qp(1) = .1/3
C      qp(2) = .01
C      qp(3) = .02
C      qp(1) = .125d0/1
C      qp(2) = .125d0*0
C      qp(3) = .125d0*0

      nk(1) = n1
      nk(2) = n2
      nk(3) = n3

C ... Find four mesh points closest to qp.  qqb is proj of qp onto qb
      call dinv33(qb,0,qbi,det)
      call dgemm('N','N',3,1,3,1d0,qbi,3,qp,3,0d0,qqb,3)
C     print *,qqb
      call iinit(ii,12)
      do  10  i = 1, 3
        if (qqb(i) < 0)     qqb(i) = qqb(i) + nk(i)
        if (qqb(i) > nk(i)) qqb(i) = qqb(i) - nk(i)
        ii(4,i) = int(qqb(i))
        ii(i,i) = int(qqb(i))+1
        xl(i) = qqb(i) - ii(4,i)
        if (xl(i) >= .5d0) then
          ii(4,i) = ii(4,i)+1
          ii(i,i) = ii(4,i)-1
        endif
        do  11  m = 1, 3
   11   if (i /= m) ii(m,i) = ii(4,i)
   10 continue
C ... ii is needed to extract function values at four corners
      do  14  i = 1, 3
        qii(i) = ii(i,i) - ii(4,i)
        qpmq4(i) = qqb(i) - ii(4,i)
        do  12  m = 1, 4
   12   ii(m,i) = mod(ii(m,i),nk(i))+1
   14 continue
      do  15  m = 1, 4
   15 fm(m) = f(ii(m,1),ii(m,2),ii(m,3))

      fqp = fm(4)
      do  20  i = 1, 3
   20 fqp = fqp + qpmq4(i)*(fm(i)-fm(4))*qii(i)

C      do  98  m = 1, 4
C   98 print *, (ii(m,i), i = 1,3)
C      print *, (qii(i), i = 1,3)
C      print *, (qpmq4(i), i = 1,3)
C      print *, (sngl(fm(m)), m=1,4)
C      print *, fqp
C      stop

      end
