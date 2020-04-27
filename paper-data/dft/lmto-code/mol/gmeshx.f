      subroutine gmeshx(alat,plat,iq,ecut,n1,n2,n3,i1p1,n,np,
     .  nbp,ibp,bp,idx,p)
C- Tabulates G vectors centered around the origin or the corresponding
C  points in real space (iq=1).  Reciprocal and real-space
C  lattice vectors are scaled by alat, with Pi . Gj = 2 pi delta(i,j)
C  Points are tabulated as multiples of G1,G2,G3, using integers
C  (0,1,2,..n/2,-n/2,-n/2+1,...-1) for each n = n1,n2,n3 (n odd)
C  (0,1,2,..(n-1)/2,-n/2,-n/2+1...-1) for each n = n1,n2,n3 (n even)
C  Real space vectors are for integer multiples of P1/n1, P2/n2, P3/n3.
C  Lattice vectors stored in p;  n used only for dimensioning p.
C  np,idxq,ecut: if np eq n, these quantities are not used.  Otherwise:
C  Recip vectors G in au; those for which G*G > ecut are discarded.
C  idx(1..n) maps i=1..n to ip=1..np.  Thus p(idx(i))
C  On output, np is the number of points not discarded
C  Array ibp points to members of p that are ;boundary points; bp
C  holds vector to be added to p to put G at the opposite boundary
      implicit none
      integer n,n1,n2,n3,i1p1,idx(1),iq,np,nbp
      double precision alat,plat(3,3),p(n,*),bp(n,*)
      double precision glat(3,3),vol,ecut
      integer i1,i2,i3,j1,j2,j3,m1,m2,m3,i,ibp(*),iprint,k1,k2,k3
      logical lbp

C --- Setup ---
      i1 = i1p1-1
      if (iq == 1) then
        do  10  i = 1, 3
          glat(i,1) = alat*plat(i,1)/n1
          glat(i,2) = alat*plat(i,2)/n2
          glat(i,3) = alat*plat(i,3)/n3
   10   continue
      else
        call dinv33(plat,2,glat,vol)
        do  20  i = 1, 9
   20   glat(i,1) = glat(i,1)/alat
      endif

C --- For each mesh point, do ---
      m1 = (n1+1)/2
      m2 = (n2+1)/2
      m3 = (n3+1)/2
      i = 0
      np = 1
      nbp = 1
      do  30  i3 = 0, n3-1
        j3 = i3 - n3*(i3/m3)
        do  30  i2 = 0, n2-1
        j2 = i2 - n2*(i2/m2)
c       do  30  i1 = 0, n1-1
        j1 = i1 - n1*(i1/m1)
        i = i+1
        p(np,1) = j1*glat(1,1) + j2*glat(1,2) + j3*glat(1,3)
        p(np,2) = j1*glat(2,1) + j2*glat(2,2) + j3*glat(2,3)
        p(np,3) = j1*glat(3,1) + j2*glat(3,2) + j3*glat(3,3)

C ---   Check for boundary point ---
        ibp(nbp) = np
        k1 = 0
        k2 = 0
        k3 = 0
        if (n1+j1+j1 == 0) k1 = n1
        if (n2+j2+j2 == 0) k2 = n2
        if (n3+j3+j3 == 0) k3 = n3
        lbp = .false.
        if (k1+k2+k3 /= 0) then
          bp(nbp,1)= k1*glat(1,1)+ k2*glat(1,2)+ k3*glat(1,3)+ p(np,1)
          bp(nbp,2)= k1*glat(2,1)+ k2*glat(2,2)+ k3*glat(2,3)+ p(np,2)
          bp(nbp,3)= k1*glat(3,1)+ k2*glat(3,2)+ k3*glat(3,3)+ p(np,3)
          if (bp(nbp,1)**2+bp(nbp,2)**2+bp(nbp,3)**2 <= ecut) then
C            print 345, nbp, bp(nbp,1), bp(nbp,2), bp(nbp,3),
C     .        bp(nbp,1)**2+bp(nbp,2)**2+bp(nbp,3)**2
C  345       format(' gmeshx: bound pnt',i4,4f11.6,6i3)
            nbp = nbp+1
            lbp = .true.
          endif
        endif
C ---   Accept point ---
        if (p(np,1)**2+p(np,2)**2+p(np,3)**2 <= ecut) then
C          print 333, i, p(np,1), p(np,2), p(np,3),
C     .      p(np,1)**2+p(np,2)**2+p(np,3)**2,i1,i2,i3,j1,j2,j3,k1,k2,k3
C  333     format(' gmeshx:  take pnt',i4,4f11.6,6i3,'  :',3i3)
          idx(i) = np
          np = np+1
C ---   Reject point ---
        else
          idx(i) = 0
C          print 334, i, p(np,1), p(np,2), p(np,3),
C     .      p(np,1)**2+p(np,2)**2+p(np,3)**2,i1,i2,i3,j1,j2,j3
C  334     format(' gmeshx: chuck pnt',i4,4f11.6,6i3)
C ---   If original point rejected, but boundary point not, swap ---
          if (lbp) then
            nbp = nbp - 1
            p(np,1) = bp(nbp,1)
            p(np,2) = bp(nbp,2)
            p(np,3) = bp(nbp,3)
            idx(i) = np
C           print *, ' ... swap bound point',np,' into mesh point'
            np = np+1
          endif
        endif
   30 continue
      np = np-1
      nbp= nbp-1

C --- Append boundary points to p ---
      do  40  i = 1, 3
   40 call dpcopy(bp(1,i),p(np+1,i),1,nbp,1d0)

      if (iprint() >= 60) print 335, np, nbp
  335 format('gmeshx: np=',i4,'  nbp=',i4)

      end
