      subroutine pgfg2g(ld1,ld2,ldiag,offpi,offpj,dpf,ddpf,gij,ierr)
C- Convert g_ij = (P-S)^-1 to G_ij by energy scaling
C-----------------------------------------------------------------------
Ci Inputs
Ci   ld1,ld2: dimension of g_ij
Ci   ldiag:   T add diagonal term -1/2 dotdot-P/dot-P
Ci   offpi,offpj: offsets to dpf and ddpf
C-----------------------------------------------------------------------
      implicit none
      integer offpi,offpj,ld1,ld2,ierr
      logical ldiag
      double precision gij(ld1,ld2,2),ddpf(2,*)
      double complex dpf(1),xxc
      integer i,j

c      ld0 = pgplp(4,max(min(ipl,npl-1),0))
c      offpi = isum(ipl,pgplp(4,0),6)

C      call zprm('dotP',2,dpf(1+offpi),ld1,ld1,1)
C      call zprm('dotdotP',2,ddpf(1,1+offpi),ld1,ld1,1)

C      if (ldiag) then
C        ierr = 0
C        do  8  i = 1, ld1
C          if (gij(i,i,2) < 0) ierr = ierr+1
C          if (gij(i,i,2) < 0) print *, 'g',i, gij(i,i,2)
C   8   continue
C      endif

C --- g_ij <- sqrt(P^dot_i) g_ij sqrt(P^dot_j) ---
      do  10  j = 1, ld2
      do  10  i = 1, ld1
        xxc = sqrt(dpf(i+offpi)*dpf(j+offpj))
     .       *dcmplx(gij(i,j,1),gij(i,j,2))
        gij(i,j,1) = dble(xxc)
        gij(i,j,2) = dimag(xxc)
   10 continue

C --- g_ij <- g_ij - 1/2 P^dot_i / P^dotdot_i ---
      if (ldiag) then
        do  20  i = 1, ld1
          gij(i,i,1) = gij(i,i,1) + ddpf(1,i+offpi)
          gij(i,i,2) = gij(i,i,2) + ddpf(2,i+offpi)
   20   continue
      endif

C ... debugging
      if (ldiag) then
        ierr = 0
        do  30  i = 1, ld1
          if (gij(i,i,2) > 0) ierr = ierr+1
          if (gij(i,i,2) > 0) print *, i, gij(i,i,2)
   30   continue
      endif

      end
