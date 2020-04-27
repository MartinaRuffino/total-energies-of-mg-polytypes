      subroutine srvsym(ngrp,g,ag,tau,plat,nR,Rv,Rmax,ips0,bgv)
C- Setup for symmetrization of a function on a real-space lattice
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ngrp,g,ag   space group
Ci   nR,Rv       lattice vectors, sorted by len(Rv-tau)
Ci   tau         lattice of points on tau+Rv
Ci   plat        primitive lattice vectors, in units of alat
Ci   nR          number of lattice translation vectors
Ci   Rv          lattice translation vectors
Ci   Rmax        If Rmax>0, ignore rotated vectors with length > Rmax
Ci               Also, do not assume vectors sorted by length.
Co Outputs:
Co   ips0        pointer to first vector in star of this vector
Co   bgv         phase factor sum; see Remarks
Cr Remarks
Cr   Prescription for finding equivalent vectors
Cr     Each vector (Rv-tau) is rotated by space group operation to
Cr     g (Rv-tau) + ag.  The following criteria must be satisfied
Cr     for the rotated vector to be included as equivalent to the first:
Cr       1.  The translation of the origin must be a lattice vector.
Cr           (i.e. ag must be a multiple of a lattice vector)
Cr           This excludes vectors whose origin starts at some
Cr           different element of the same class.
Cr       2.  T' = g (Rv-tau) + ag + tau must be a lattice vector.
Cr           This again excludes vectors connecting to different
Cr           different element of the same class.
Cr       3.  (this should no longer be needed; keep for safety)
Cr           When tau is nonzero, the lattice translation vectors that
Cr           define the rotated tau in its symmetry- related positions
Cr           may have different lengths, some of which may be larger
Cr           than those in the list of vectors Rv.  Input Rmax>0 sets
Cr           the length of the longest vector in list.
Cr           If Rmax<0, these vectors are flagged by ips0(:)<0
Cu Updates
Cu   24 Jun 09 Flag vectors that do not have a complete star
Cu   6 Sep 00 fixed some bugs
C ----------------------------------------------------------------------
      implicit none
      integer ngrp,nR,ips0(nR)
      double precision g(3,3,ngrp),ag(3,ngrp),Rv(nR,3),plat(3,3),tau(3),
     .  Rmax
      double complex bgv(nR)
      integer i,i00,irep,i0,nstar,k,j,j0,iprint,ksum,kstar,m,lgunit,
     .  nskip
      double precision df,gg0,gg,fac,vv,qlat(3,3),
     .  v(3),v2(3),xx,rsiz,bsum,vsize,rmax2,fuzz
      parameter (fuzz=1d-8)
      logical latvec,lok
      rsiz(m) = (Rv(m,1)-tau(1))**2 +
     .          (Rv(m,2)-tau(2))**2 +
     .          (Rv(m,3)-tau(3))**2

C     tpi = 8d0*datan(1d0)
      call dinv33(plat,1,qlat,xx)
      rmax2 = rmax*rmax

C ... Setup: initialize ips0
      do  10  i = 1, nR
        ips0(i) = 0
        bgv(i) = (0d0,0d0)
   10 continue

C --- Main loop: look for next unclassified vector ---
      i00 = 1
      do  20  irep = 1, nR+1
        i0 = 0
        do  22  i = i00, nR
          i0 = i
          if (ips0(i) == 0) goto 80
   22   continue
        goto 81
   80   continue

C   --- Apply all point ops, find in list, add to phase sum ---
        nstar = irep
        do  30  k = 1, ngrp

C     ... T'-tau = g (T-tau) + ag => R' = g R + ag + tau
          do  32  i = 1, 3
   32     v(i) = g(i,1,k)*(Rv(i0,1)-tau(1)) +
     .           g(i,2,k)*(Rv(i0,2)-tau(2)) +
     .           g(i,3,k)*(Rv(i0,3)-tau(3)) + ag(i,k) + tau(i)

C         Exclude sites not connected to tau by lattice vector, ie
C         vectors connecting R to a different member of the class.
C         Such vectors are not among the class of vectors tau+T
C         available to symmetrize
          lok = latvec(1,fuzz,qlat,v)
C         Exclude shifts of origin
          lok = lok .and. latvec(1,fuzz,qlat,ag(1,k))

C     ... Find index j0 to lattice vector corresponding to this v
          if (lok) then

            do  35  j = i0, nR
              v2(1) = v(1) - Rv(j,1)
              v2(2) = v(2) - Rv(j,2)
              v2(3) = v(3) - Rv(j,3)
              df = v2(1)**2 + v2(2)**2 + v2(3)**2
              j0 = j
              if (df < fuzz) goto 70
   35       continue


C       ... No lattice vector found
            if (ips0(i0) < 0) goto 30
            if (rmax < 0) then
            do  i = i0+1, nR
              if (ips0(i) == i0) ips0(i) = -i0
            enddo
            ips0(i0) = -ips0(i0)
            goto 30
            endif

C           No lattice vector found.  No error if R' > Rmax.
C           NB: there still may be a problem because some Rv at
C           boundaries may be included while others not.
            if (rmax > 0) then
              vsize = (v(1)-tau(1))**2+(v(2)-tau(2))**2+(v(3)-tau(3))**2
              if (vsize > rmax2-fuzz) goto 30
            endif

            print 601, i0,k,Rv(i0,1),Rv(i0,2),Rv(i0,3),v
  601       format(' ---- vec',i6,'   op',i3,2x,3f8.4,'  ->',3f8.4)
            call rxi('SRVSYM: cannot find mapped vector in list:',i0)
   70       continue
            if (ips0(i0) < 0) then
              ips0(j0) = -i0
            else
              ips0(j0) = i0
            endif

            if (ips0(j0) > 0) then
C             scalp = Rv(j0,1)*ag(1,k)+Rv(j0,2)*ag(2,k)+Rv(j0,3)*ag(3,k)
C             bgv(j0) = bgv(j0) + cdexp(dcmplx(0d0,tpi*scalp))
              bgv(j0) = bgv(j0) + 1

C             call daxpy(3,-1d0,tau,1,v,1)
              if (iprint() > 90) print 603, i0,j0,k,
     .          Rv(i0,1)-tau(1)*0,Rv(i0,2)-tau(1)*0,Rv(i0,3)-tau(1)*0,v
  603       format(' vec',i5,'  ->',i5,'   op',i3,2x,3f8.4,'  ->',3f8.4)

            else

              if (iprint() > 90)
     .          print 602, i0,Rv(i0,1),Rv(i0,2),Rv(i0,3),k
  602         format(' vec',i6,3f8.4,'   op',i3,' missing star')

            endif

          endif
   30   continue
        i00 = i0
   20 continue
      call rxi('SRVSYM: this cannot happen. irep=',irep)
   81 continue

C      if (iprint() >= 20) print 878, nR,nstar,tau
C  878 format(/' Srvsym:',i6,'  vectors in',i6,' stars.  tau =',3f12.6)
      if (iprint() >= 30) call awrit3(
     .  '%N srvsym:  %i vectors in %i stars.  tau =%3:1,3;3d',' ',80,
     .  lgunit(1),nR,nstar,tau)

C --- Multiply phase sums by (star order)/(eff group order) ---
      nskip = 0
      ksum = 0
      do  i0 = 1, nR
        if (ips0(i0) < 0) then
          nskip = nskip + 1
        endif
        if (ips0(i0) == i0) then
          kstar = 0
          bsum = 0d0
          gg0 = rsiz(i0)
          do  i = i0, nR
            if (ips0(i) == i0) then
              kstar = kstar+1
              bsum = bsum + dble(bgv(i))
            endif
            if (rmax == 0) then
              gg = rsiz(i)
              if (dabs(gg-gg0) > 1d-3) goto 78
            endif
          enddo
   78     continue
          ksum = ksum+kstar
          fac = dble(kstar)/bsum
          do  i = i0, nR
            if (ips0(i) == i0) bgv(i) = bgv(i)*fac
            if (rmax == 0) then
              gg = rsiz(i)
              if (dabs(gg-gg0) > 1d-3) goto 79
            endif
          enddo
   79     continue
        endif
      enddo
      if (ksum+nskip /= nR) call rxi('SRVSYM error, ksum=',ksum)

C --- Printout ---
      if (iprint() < 50) return
      j = min(nR,300)
      if (iprint() >= 60) j = nR
      nstar = 0
      do  50  i0 = 1, j
        if (ips0(i0) == i0) then
          nstar = nstar+1
          vv = dsqrt(rsiz(i0))
          print 250, nstar,vv
  250     format(/' Star',i5,'     vector length',f10.4)
          kstar = 0
          do  56  i = i0, nR
            if (ips0(i) == i0) then
              kstar = kstar+1
              print 251, kstar,i,Rv(i,1),Rv(i,2),Rv(i,3),bgv(i)
  251         format(i5,i8,2x,3f7.2,3x,2f7.2)
            endif
   56     continue
        endif
   50 continue

      end
      subroutine shfspc(g,ag,nsgrp,tau,ib,plat,agT)
C- Shift translation part of space group corresponding to lattice shift
C ----------------------------------------------------------------------
Ci Inputs
Ci Inputs
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   nsgrp :number of space group operations
Ci   tau,ib:lattice translation: pos(i) translated to pos-tau(ib)
Co Outputs
Ci   agT   :translation part of space group for translated lattice
Cr Remarks
Cr   For group op (g,ag) centered at 0, the corresponding group op
Cr   centered at tau is (gT,agT) = (g,ag+g*tau-tau)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer nsgrp,ib
      double precision g(3,3,nsgrp),ag(3,nsgrp),agT(3,nsgrp),tau(3,ib)
      double precision plat(3,3)
      integer i,k
      double precision v(3),qlat(3,3),xx

      call dinv33(plat,1,qlat,xx)

C ... Make agT = a + g tau - tau
      do  10  k = 1, nsgrp
      call dmpy(g(1,1,k),3,1,tau(1,ib),3,1,v,3,1,3,1,3)
      do  20  i = 1, 3
   20 agT(i,k) = v(i) + ag(i,k) - tau(i,ib)
      call shorbz(agT(1,k),agT(1,k),plat,qlat)
   10 continue
C     call prmx('agT',agT,3,3,nsgrp)
      end
C      subroutine snot(g,ag,ng,pos,plat)
C      implicit none
C      integer ng
C      double precision g(3,3,ng),ag(3,ng),pos(3),plat(3,3)
C      double precision v2(3),pos2(3,ng),ddot,qlat(3,3),v(3),xx
C
C      double precision pos3(3,ng)
C      integer i,k
C
C      print 333, pos
C  333 format('snot: pos=',3f12.6)
C      call dinv33(plat,1,qlat,xx)
C
C      do  10  k = 1, ng
C        call dmpy(g(1,1,k),3,1,pos,3,1,v2,3,1,3,1,3)
C        do  31  i = 1, 3
C          pos2(i,k) = v2(i) + ag(i,k)
C          pos3(i,k) = v2(i) + ag(i,k) - pos(i)
C   31   continue
C   10 continue
C
C      do  20  k = 1, ng
C        do  38  i = 1, 3
C        v(i)  = qlat(1,i)*pos3(1,k)+qlat(2,i)*pos3(2,k)+
C     .          qlat(3,i)*pos3(3,k)
C        v2(i) = qlat(1,i)*pos2(1,k)+qlat(2,i)*pos2(2,k)+
C     .          qlat(3,i)*pos2(3,k)
C   38 continue
C
CC        print 603, (pos2(i,k),i=1,3),
CC     .    ddot(3,pos2(1,k),1,pos2(1,k),1),(pos3(i,k),i=1,3)
C        print 603, k,(pos2(i,k),i=1,3),
C     .    ddot(3,pos2(1,k),1,pos2(1,k),1), v,v2
C  603   format(i3,3f12.8,1x,f10.4,3x,3f12.8,3x,3f12.8)
C   20 continue
C      stop
C
C      end
