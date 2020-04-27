      subroutine lgen(pqlat,bmax,nv,nvmax,vecs)
C- Generate lattice vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   pqlat:primitive lattice vectors
Ci   bmax: largest radius for vector?
Ci   nvmax: maximum number of vectors allowed
Co Outputs
Co   nv,vecs
Cr Remarks
Cr   Lattice vectors are sorted by increasing length.
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nvmax
      double precision pqlat(3,3),vecs(3,0:*),bmax
C Local parameters
      double precision bmax2,v2
      integer i,j,k,imax,jmax,kmax,nv,m,ivck(3),iprint

      ivck(1) = 1
      ivck(2) = 1
      ivck(3) = 1
      call latlim(pqlat,bmax,imax,jmax,kmax)
      imax = max(imax,1)
      jmax = max(jmax,1)
      kmax = max(kmax,1)
      bmax2 = bmax*bmax
      nv = 0
      do  i = -imax, imax
      do  j = -jmax, jmax
      do  k = -kmax, kmax
        v2 = 0
            do  m = 1, 3
          vecs(m,nv) = i*pqlat(m,1) + j*pqlat(m,2) + k*pqlat(m,3)
          v2 = v2 + vecs(m,nv)**2
            enddo
C       --- Include as lattice vector if a plat or length > bmax**2 ---
        if (iabs(i) + iabs(j) + iabs(k) == 1 .and. v2 <= bmax2) then
          if (i == 1) ivck(1) = 0
          if (j == 1) ivck(2) = 0
          if (k == 1) ivck(3) = 0
        endif
C   --- Case a lattice vector found ---
        if (v2 <= bmax2) then
          nv = nv+1
          if (nv >= nvmax) call rx(' lgen: too many vectors')
        endif
      enddo
      enddo
      enddo
      if (ivck(1)+ivck(2)+ivck(3) /= 0 .and. iprint() >= 20) then
        print 333, ivck
  333   format(/' lgen: added missing plat: ivck=',3i2)
        if (3*nv > nvmax)call fexit(-1,9,' lgen: too many vectors',0)
        if (ivck(1)+ivck(2)+ivck(3) /= 1)
     .    call fexit(-1,9,' lgen: more than 1 missing plat',0)
        do  m = 1, 3
          v2 = ivck(1)*pqlat(m,1)+ivck(2)*pqlat(m,2)+ivck(3)*pqlat(m,3)
          call dcopy(nv,vecs(m,1),3,vecs(m,nv+1),3)
          call dcopy(nv,vecs(m,1),3,vecs(m,2*nv+1),3)
          call daxpy(nv, 1d0,v2,0,vecs(m,nv+1),3)
          call daxpy(nv,-1d0,v2,0,vecs(m,2*nv+1),3)
        enddo
        nv = 3*nv
      endif
      call bsortl(vecs,nv)
      end
