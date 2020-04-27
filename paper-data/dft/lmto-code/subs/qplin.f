      integer function qplin(qpi,plat,nqfbz,qfbz,ipq,iqfbz,wt)
C- Find four QP nearest to given up, and linear interpolation weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   qpi   :qp for which to find interpolation weights
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nqfbz :number of qp in full BZ
Ci   qfbz  :qp in full BZ
Ci   ipq   :ipq(1)=0 => not used
Ci         :ipq(i1,i2,i3) points to the irreducible qp into which
Ci         :mesh point (i1,i2,i3) is mapped (bzmesh.f)
Co Outputs
Co   iqfbz :iqfbz(:,1): indices to qp in BZ to be used in interpolation
Co         :iqfbz(:,2): not used if ipq(1)=0
Co         :            else, points in irr of BZ equivalent to iqfbz(:,1)
Co   wt    :interpolation weights
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr  Linear interpolation with barycentric coordinates for tetrahedra
C   Linearly Interpolate from pts at corners of tetrahedra
Cu Updates
Cu   25 Jul 12 Improved handling of singular cases
Cu   18 Jul 12 First created
C ----------------------------------------------------------------------
      implicit none
      integer nqfbz,ipq(*),iqfbz(4,2)
      double precision qpi(3),qfbz(3,nqfbz),plat(3,3),wt(4)
C ... Local parameters
      integer,allocatable :: iprm(:)
      real(8),allocatable:: qp(:,:)
      integer i,j,nsim,ipvt(4),iq4
      double precision swt,qlat(3,3),xx,det(2),wk(16),ddot
      double precision norm(4,4),norm0(4,4)
      real(8):: tolq=1d-7,ltolq=-7d0
      procedure(real(8)) :: angle, dlength

      if (nqfbz < 4) goto 99

      call mkqlat(plat,qlat,xx)

C     call shorps(1,qlat,(/61,61,61/),qpi,q)

C ... Sort (qfbz-qpi) to get 4 k-points in BZ nearest to qp
      allocate(qp(3,nqfbz),iprm(nqfbz))
      do  i = 1, nqfbz
        call shorbz(qfbz(:,i)-qpi(:),qp(1,i),qlat,plat)
      enddo
      call dvheap(3,nqfbz,qp,iprm,tolq,11)

C ... If first point coincides with a mesh point, wt(1)=1
      if (dlength(3,qp(1,iprm(1)),1) < tolq) then
        norm(1,:) = 0
        norm(1,1) = 1
        nsim = 4
        goto 20
      endif

C ... If points 1,2,3 are collinear, find point that is not
      iq4 = 4
      do  while (abs(abs(angle(3,qp(:,iprm(1))-qp(:,iprm(2)),
     .               qp(:,iprm(1))-qp(:,iprm(3))))-1) < tolq)
        iq4 = iq4+1
        if (iq4 > nqfbz) goto 99
        i = iprm(3)
        iprm(3) = iprm(iq4)
        iprm(iq4) = i
      enddo

C ... Find weights for each of four points
C     Make nsim = # independent pts; wt = weights for points
C     Weights generated in norm(1,1..nsim)
      nsim = 4
      wt = 0
C     norm(1:4,2:4) <- transpose of qp
      call dmcpy(qp,1,3,norm(1,2),4,1,4,3)
      do  i = 1, 4
        norm(i,2:4) = qp(1:3,iprm(i))
      enddo
C     norm(1:4,1) <- 1
      call dvset(norm,1,4,1d0)
      call dcopy(16,norm,1,norm0,1)
      call icopy(4,iprm,1,iqfbz,1)
C     call prmx('norm',norm,4,nsim,nsim)
      goto 12

C ... Shuffle 4th point, in case |det| is too small
   10 continue
      iq4 = iq4+1
      if (iq4 > nqfbz) goto 99
      call dcopy(16,norm0,1,norm,1)
      i = iprm(4)
      iprm(4) = iprm(iq4)
      iprm(iq4) = i
      norm(4,2:4) = qp(1:3,iprm(4))
C     call prmx('norm',norm,4,nsim,nsim)

   12 continue
      call dgefa(norm,4,nsim,ipvt,j)
      if (j /= 0) goto 10
      call dgedi(norm,4,nsim,ipvt,det,wk,11)
      if (det(2) < ltolq) goto 10
C     call prmx('norm',norm,4,nsim,nsim)
C     If inversion failed, exit
      if (j /= 0) goto 99

   20 continue
C     If any weight exceeds 10 or swt<>1, exit
      swt = 0
      do  j = 1, nsim
        iqfbz(j,1) = iprm(j)
        if (ipq(1) /= 0) then
          iqfbz(j,2) = ipq(iprm(j))
        endif
        wt(j) = norm(1,j)
        if (wt(j) > 10) cycle
        swt = swt + wt(j)
      enddo
      if (abs(swt-1) > tolq) goto 99
C     All tests passed ... use this set of points
C     enddo
C     call prmx('inverse',norm,4,nsim,nsim)
      if (abs(swt-1) > tolq) goto 99

C ... Printout header information
      call info5(40,0,0,' Interp SE for qpi=%3:1,5;5d, '//
     .  'wt=%n:1,5;5d',qpi,nsim,wt,0,0)
      do  j = 1, nsim
        xx = dsqrt(ddot(3,qp(1,iprm(j)),1,qp(1,iprm(j)),1))
        call daxpy(3,1d0,qpi,1,qp(1,iprm(j)),1)
        if (ipq(1) /= 0) then
          call info8(41,0,0,' qi(%i) =%3;12,6D  |qi-qpi|=%,6;6d'//
     .      '  iq=%i irr=%i  wt=%d',j,qp(1,iprm(j)),xx,iqfbz(j,1),iqfbz(j,2),
     .      wt(j),0,0)
        else
          call info5(41,0,0,' qi(%i) =%3;12,6D  |qi-qpi|=%,6;6d'//
     .      '  iq=%i wt=%d',j,qp(1,iprm(j)),xx,iqfbz(j,1),wt(j))
        endif
      enddo
      qplin = 0
      return

C ... Error exit
   99 continue
      call info2(20,0,0,' wt gen failed for qp=%3:1,5;5d, '//
     .  'wt =%4:1;9F',qpi,wt)
      qplin = -1

      end

      double precision function angle(n,a,b)
      implicit none
      integer n
      double precision a(n),b(n)
      procedure(real(8)) :: ddot

      angle = ddot(n,a,1,b,1)/dsqrt(ddot(n,a,1,a,1)*ddot(n,b,1,b,1))
      end
