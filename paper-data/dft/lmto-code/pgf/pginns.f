      subroutine pginns(ld,n,r,rvec,dlnr,wt,lg00,gi00,igb00,sii,mode,
     .  wk,lgnn,ginn)
C- Make gs_nn given gs_00 and semi-infinite BC on rhs
C ----------------------------------------------------------------
Ci  Inputs:
Ci    n    : make gi_n,n (see mode)
Ci    igb00: inverse of bulk GF, i.e. (g^B_00)^-1 (mode=0,1,3)
Ci    gi00 : interface GF at 0
Ci    ginn : interface GF at n (mode=2)
Ci    dlnr :
Ci    wt   : used in adiabatic approximation
Ci    mode : ones digit
Ci         : 0
Ci         : 1 put info into rvec needed for backward updates
Ci         : 2 Use info from rvec to make backward updates
Ci         :   (see mode, Outputs)
Ci         : 3 generate information needed for g0n,gn0
Ci    wk   : double complex work array of dimension ld*ld
Ci    ginn : (mode 2)
Cio Inputs/Outputs
Cio   rvec :  parameters related to eigenvectors of normal modes
Cio        : mode 1 input: eigenvectors
Cio        :   (1,1)  Q     (1,2) P
Cio        :   (2,1): Q^-1  (2,2) P^-1
Cio        : rvec is OVERWRITTEN on mode 1 output
Cio        : mode 1 output and mode 2 input:
Cio        :   (1,1) 1+A      (2,1) gi_0n gi_nn^-1
Cio        :   (1,2) gi_n0    (2,2) gi_00 (g^B_00)^-1
Cio        : mode 3 output:
Cio        :   (1,1) 1+A               (2,1) gs_0n
Cio        :   (1,2) P r^n P^-1 gi00   (2,2) gs_n0
Cio   sii  : On input
Cio        : sii(2) QrQ^-1   sii(3) PrP^-1 (see pgbulk)
Cio        : On output:
Cio        : sii(1) C = (B-1) (1 + A B)^-1  CHECK
Co  Outputs
Co   ginn  :(mode=0,1,3) surface diagonal GF at layer n
Co         :(mode=2)     proper  diagonal GF at layer n
Co   gi00  :(mode=2)     proper  diagonal GF at layer 0
Cr  Remarks
Cr    Let A = P r^n P^-1 (g^i_00 g^B_00^-1 - 1) Q r^n Q^-1
Cr    Then
Cr      g^i_n,n     = (1 + A) g^B_00
Cr      g^i_n+1,n   = PrP^-1 (1 + A) g^B_00 = PrP^-1 g^i_n,n
Cr      g^i_n,n+1   = (1 + A) QrQ^-1 g^B_00
Cr      g^i_n+1,n+1 = (1 + PrP^-1 A QrQ^-1 ) g^B_00
Cr    Also let B = QrQ^-1 PrP^-1  and C = (B-1) (1 + A B)^-1
Cr    g^i_nn = (1+A) g^B_00
Cr    g^s_nn = [1 - (1+A) (B^-1 + A)^-1 ] (1+A) g^B_00
Cr           = (1-B) (1 + AB)^-1 (1+A) g^B_00
Cr           = (1-B) (1 + AB)^-1 g^i_n,n = -C g^i_n,n
Cr    g^i_0n = g^i_00 (g^B_00)^-1 Q r^n Q^-1 (g^B_00)
Cr    g^s_0n = g^i_0n (g^i_nn)^-1 g^s_nn
Cr           = g^i_00 (g^B_00)^-1 Q r^n Q^-1 (1+A)^-1 g^s_nn
Cr    g^i_n0 = P r^n P^-1 g^i_00
Cr    g^s_n0 = -C g^i_n,0 = -C P ^n P^-1 gi00
C ----------------------------------------------------------------
      implicit none
      integer ld,lg00,n,lgnn,mode
      double precision igb00(ld,ld,2),gi00(ld,ld,2),ginn(ld,ld,2),
     .  wk(ld,ld,2),rvec(ld,2,ld,2,2),r(ld,2,2),dlnr(ld,2,2),wt(n),
     .  sii(ld,ld,3,2)
      integer i,j,ierr,mode0,mode1
C     parameter (PRTG=80)
      logical lback

C     call wkchk('test')
      call rxx(mod(lgnn,2) == 0,'pginns: no memory allocated for ginn')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      if (mode0 == 2) goto 100
      lback = mode0 == 1 .or. mode0 == 3

C --- sii(1) <- B = QrQ^-1 PrP^-1 ---
      call yygemm('N','N',ld,ld,ld,1d0,sii(1,1,2,1),sii(1,1,2,2),ld,
     .  sii(1,1,3,1),sii(1,1,3,2),ld,0d0,sii(1,1,1,1),sii(1,1,1,2),ld)
C     call yprm('B',2,sii,ld*ld*3,ld,ld,ld)
C ... Now sii(2..3) free for workspace

C --- Make A as defined in Remarks; put in rvec(1,1) ---
C ... sii3 <- P r^n P^-1
      call pgprp(n,mode1*10+2,ld,r,rvec,dlnr,wt,wk,ld,sii(1,1,3,1),ld*3)
C ... wk <- (g00^i g^B_00^-1 - 1) ...
      call yygemm('N','N',ld,ld,ld,1d0,gi00,gi00(1,1,2),ld,
     .  igb00,igb00(1,1,2),ld,0d0,wk,wk(1,1,2),ld)
C ... (lback): rvec(2,2) <- g00^i g^B_00^-1, rvec(1,2) <- Pr^nP^-1 gi00
      if (lback) then
        do  12  j = 1, ld
        do  12  i = 1, ld
        rvec(i,2,j,2,1) = wk(i,j,1)
   12   rvec(i,2,j,2,2) = wk(i,j,2)
        call yygemm('N','N',ld,ld,ld,1d0,sii(1,1,3,1),sii(1,1,3,2),ld,
     .    gi00,gi00(1,1,2),ld,0d0,rvec(1,1,1,2,1),rvec(1,1,1,2,2),ld*2)
      endif
      do  10  i = 1, ld
   10 wk(i,i,1) = wk(i,i,1) - 1
C ... sii2 <- P r^n P^-1 (g00^i g^B_00^-1 - 1) ...
      call yygemm('N','N',ld,ld,ld,1d0,sii(1,1,3,1),sii(1,1,3,2),ld,
     .  wk(1,1,1),wk(1,1,2),ld,0d0,sii(1,1,2,1),sii(1,1,2,2),ld)
C ... sii3 <- Q r^n Q^-1
      call pgprp(n,mode1*10+1,ld,r,rvec,dlnr,wt,wk,ld,sii(1,1,3,1),ld*3)
C ... (lback): rvec(2,1) = g00^i g^B_00^-1 Q r^n Q^-1 ...
      if (lback) then
        call yygemm('N','N',ld,ld,ld,1d0,rvec(1,2,1,2,1),
     .    rvec(1,2,1,2,2),
     .    ld*2,sii(1,1,3,1),sii(1,1,3,2),ld,0d0,
     .    rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2)
      endif
C ... rvec(1,1) <- A = P r^i P^-1 (g00^i g^B_00^-1 - 1) Q r^i Q^-1 ...
      call yygemm('N','N',ld,ld,ld,1d0,sii(1,1,2,1),sii(1,1,2,2),ld,
     .  sii(1,1,3,1),sii(1,1,3,2),ld,0d0,rvec,rvec(1,1,1,1,2),ld*2)
C     call yprm('A',2,rvec,(ld*2)**2,ld*2,ld,ld)

C --- Make (B-1) (1 + A B)^-1 ---
C ... sii(2) <- 1+AB, sii(1) <- B-1, sii(1) <- (B-1) (1 + A B)^-1
      call yygemm('N','N',ld,ld,ld,1d0,rvec,rvec(1,1,1,1,2),ld*2,
     .  sii(1,1,1,1),sii(1,1,1,2),ld,0d0,sii(1,1,2,1),sii(1,1,2,2),ld)
      do  20  i = 1, ld
      sii(i,i,1,1) = sii(i,i,1,1) - 1
   20 sii(i,i,2,1) = sii(i,i,2,1) + 1
C     call yprm('1+AB',2,sii(1,1,2,1),ld*ld*3,ld,ld,ld)
C ... sii(1) <- (B-1)(1+AB)-1 ...
      call yyqnvb('t',sii(1,1,2,1),sii(1,1,2,2),ld,ld,ld,
     .  wk,ld,wk,sii(1,1,1,1),sii(1,1,1,2),ld,ierr)
C     call yprm('(B-1)(1+AB)^-1',2,sii,ld*ld*3,ld,ld,ld)

C --- rvec(1,1) <- (1+A); (lback): rvec(2,1) *= (1+A)^-1 ---
      do  30  i = 1, ld
   30 rvec(i,1,i,1,1) = rvec(i,1,i,1,1) + 1
      if (lback) then
        do  32  j = 1, ld
        do  32  i = 1, ld
        sii(i,j,2,1) = rvec(i,1,j,1,1)
   32   sii(i,j,2,2) = rvec(i,1,j,1,2)
        call yyqnvb('t',sii(1,1,2,1),sii(1,1,2,2),ld,ld,ld,
     .    wk,ld,wk,rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2,ierr)
      endif

C --- rvec(1,1) <- g^i_nn = (1+A) gb00r ---
      call yyqnvb('t',igb00,igb00(1,1,2),ld,ld,ld,
     .  wk,ld,wk,rvec,rvec(1,1,1,1,2),ld*2,ierr)
C     call yprm('g^i_nn',2,rvec,(ld*2)**2,ld*2,ld,ld)

C --- Make g^s_nn = -(B-1) (1 + A B)^-1 g^i_nn = -C g^i_nn ---
      call yygemm('N','N',ld,ld,ld,-1d0,sii(1,1,1,1),sii(1,1,1,2),ld,
     .  rvec,rvec(1,1,1,1,2),ld*2,0d0,ginn,ginn(1,1,2),ld)
C     call yprm('g^s_nn',2,ginn,ld*ld,ld,ld,ld)

C --- Make off-diagonal surface g's ---
      if (mode0 == 3) then
C   ... g^s_0n = g^i_00 (g^B_00)^-1 Q r^n Q^-1 (1+A)^-1 g^s_nn
        call yygemm('N','N',ld,ld,ld,1d0,rvec(1,2,1,1,1),
     .    rvec(1,2,1,1,2),
     .    ld*2,ginn,ginn(1,1,2),ld,0d0,sii(1,1,2,1),sii(1,1,2,2),ld)
        do  34  j = 1, ld
        do  34  i = 1, ld
        rvec(i,2,j,1,1) = sii(i,j,2,1)
   34   rvec(i,2,j,1,2) = sii(i,j,2,2)
C       call yprm('g^s_0n',2,rvec(1,2,1,1,1),ld*2*ld*2,ld*2,ld,ld)
C   ... g^s_n0 = -C g^i_n0
        call yygemm('N','N',ld,ld,ld,-1d0,sii(1,1,1,1),sii(1,1,1,2),ld,
     .    rvec(1,1,1,2,1),rvec(1,1,1,2,2),ld*2,0d0,
     .    rvec(1,2,1,2,1),rvec(1,2,1,2,2),ld*2)
C       call yprm('g^s_n0',2,rvec(1,2,1,2,1),ld*2*ld*2,ld*2,ld,ld)
      endif
      call lgupac(lgnn,'p',1,1,0,7,0,0)
      call lgupac(lgnn,' pginns created',6,0,0,0,0,0)
      return

C --- Complete back-update of g00 from new ginn ---
  100 continue

C ... sii(1) <- ginn_updated (ginn_original)^-1 - 1 ..
      do  110  j = 1, ld
      do  110  i = 1, ld
      sii(i,j,1,1) = ginn(i,j,1)
      sii(i,j,1,2) = ginn(i,j,2)
      sii(i,j,2,1) = rvec(i,1,j,1,1)
  110 sii(i,j,2,2) = rvec(i,1,j,1,2)
      call yyqnvb('t',sii(1,1,2,1),sii(1,1,2,2),ld,ld,ld,
     .  wk,ld,wk,sii(1,1,1,1),sii(1,1,1,2),ld,ierr)
      do  120  i = 1, ld
  120 sii(i,i,1,1) = sii(i,i,1,1) - 1

C ... sii(2) = gi00 (gb00)^-1 Q r^n Q^-1 (1+A)^-1 sii(1) ...
      call yygemm('N','N',ld,ld,ld,1d0,rvec(1,2,1,1,1),rvec(1,2,1,1,2),
     .  ld*2,sii,sii(1,1,1,2),ld,0d0,sii(1,1,2,1),sii(1,1,2,2),ld)

C ... sii(1) = sii(2) *  (Pr^nP^-1 ginn_original)
      call yygemm('N','N',ld,ld,ld,1d0,sii(1,1,2,1),sii(1,1,2,2),ld,
     .  rvec(1,1,1,2,1),rvec(1,1,1,2,2),ld*2,0d0,sii,sii(1,1,1,2),ld)

C ... Updated gi_00 = original gi_00 - sii(1) ...
      call daxpy(ld*ld,1d0,sii,1,gi00,1)
      call daxpy(ld*ld,1d0,sii(1,1,1,2),1,gi00(1,1,2),1)
C     call yprm('updated g^i_00',2,gi00,ld*ld,ld,ld,ld)

      call lgupac(lg00,'ptmb',0,3,0,4+8,0,0)
      call lgupac(lg00,' pginns created',6,0,0,0,0,0)

      end
      subroutine pgprp(n,opt,ld,r,rvec,dlnr,wt,wk,ldw,QrQ,ldQ)
C- Make Q r^n Q^-1 or P r^n P^-1
C ----------------------------------------------------------------
Ci Inputs:
Ci   ld:   dimension of bulk PL
Ci   r:    log of the wave vectors (PRB 39, 923 (1989)).
Ci   rvec: (1,1)  Q     (1,2) P
Ci         (2,1): Q^-1  (2,2) P^-1
Ci   opt:  1 make Q r^-n Q^-1
Ci         2 make P r^n P^-1
Ci        11 make Q (r1 r2 r3 ... ) Q^-1
Co Outputs:
Co
Cr Remarks
C ----------------------------------------------------------------
      implicit none
      integer ld,opt,n,ldw,ldq
      double precision wk(ld,ldw,2),QrQ(ld,ldQ,2),rvec(ld,2,ld,2,2),
     .  r(ld,2,2),dlnr(ld,2,2),wt(n)
      double precision xxr,xxi,x2r,x2i
      integer i,j,m,mode1
      double complex lnr,dr

      m = mod(opt,10)
      mode1 = mod(opt/10,10)

C ---  wk <- Q r^n ---
C      do  10  j = 1, ld
C        xxr = r(j,m,1)
C        xxi = r(j,m,2)
C        do  12  i = 2, n
C          x2r = xxr
C          x2i = xxi
C          xxr = x2r*r(j,m,1) - x2i*r(j,m,2)
C          xxi = x2i*r(j,m,1) + x2r*r(j,m,2)
C   12   continue
C        do  14  i = 1, ld
C        wk(i,j,1) = rvec(i,1,j,m,1)*xxr - rvec(i,1,j,m,2)*xxi
C   14   wk(i,j,2) = rvec(i,1,j,m,2)*xxr + rvec(i,1,j,m,1)*xxi
C   10 continue
C
C*      call yprm('Q r^n',2,wk,ld*ldw,ld,ld,ld)
C
CC ... QrQ <- Q r^n Q^-1 ...
C      call yygemm('N','N',ld,ld,ld,1d0,wk,wk(1,1,2),ld,
C     .  rvec(1,2,1,m,1),rvec(1,2,1,m,2),ld*2,0d0,QrQ,QrQ(1,1,2),ld)
C


C ---  wk <- Q r^n ---
      do  10  j = 1, ld
        if (mode1 == 0) then
          xxr = r(j,m,1)
          xxi = r(j,m,2)
        else
          dr = dcmplx(dlnr(j,m,1),dlnr(j,m,2))
          lnr = cdexp(cdlog(dcmplx(r(j,m,1),r(j,m,2))) + wt(1)*dr)
          xxr = dble(lnr)
          xxi = dimag(lnr)
        endif
        do  12  i = 2, n
          x2r = xxr
          x2i = xxi
          if (mode1 /= 0) then
            lnr = cdexp(cdlog(dcmplx(r(j,m,1),r(j,m,2))) + wt(i)*dr)
            xxr = x2r*dble(lnr) - x2i*dimag(lnr)
            xxi = x2i*dble(lnr) + x2r*dimag(lnr)
          else
            xxr = x2r*r(j,m,1) - x2i*r(j,m,2)
            xxi = x2i*r(j,m,1) + x2r*r(j,m,2)
          endif
C         if (j == ld) print 333, i,wt(i),xxr,xxi,x2r,x2i
C  333     format(i4,f8.3,1p,4e12.4)
   12   continue
        do  14  i = 1, ld
        wk(i,j,1) = rvec(i,1,j,m,1)*xxr - rvec(i,1,j,m,2)*xxi
   14   wk(i,j,2) = rvec(i,1,j,m,2)*xxr + rvec(i,1,j,m,1)*xxi
   10 continue

C     call yprm('Q r^n',2,wk,ld*ldw,ld,ld,ld)

C ... QrQ <- Q r^n Q^-1 ...
      call yygemm('N','N',ld,ld,ld,1d0,wk,wk(1,1,2),ld,
     .  rvec(1,2,1,m,1),rvec(1,2,1,m,2),ld*2,0d0,QrQ,QrQ(1,1,2),ld)

C     call yprm('Q r^n Q^-1',2,QrQ,ld*ldQ,ld,ld,ld)

      end
