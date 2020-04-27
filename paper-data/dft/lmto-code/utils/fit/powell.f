      subroutine powell(p,xi,n,np,ftol,iter,fret)
C Passed Parameters
      implicit none
      integer iter,n,np
      double precision p(1),xi(np,1),ftol,fret
C Local Parameters
      integer i,ibig,itmax,j,nmax
      double precision del,fp,fptt,func,t
      PARAMETER (NMAX=20,ITMAX=200)
      double precision PT(NMAX),PTT(NMAX),XIT(NMAX)
      FRET=FUNC(P)
      DO 11 J=1,N
        PT(J)=P(J)
   11 CONTINUE
      ITER=0
    1 ITER=ITER+1
      FP=FRET
      IBIG=0
      DEL=0.d0
      DO 13 I=1,N
        DO 12 J=1,N
          XIT(J)=XI(J,I)
   12   CONTINUE
        CALL LINMIN(P,XIT,N,FRET)
        IF(dABS(FP-FRET).GT.DEL)THEN
          DEL=dABS(FP-FRET)
          IBIG=I
        ENDIF
   13 CONTINUE
      IF(2.d0*dABS(FP-FRET).LE.FTOL*(dABS(FP)+dABS(FRET)))RETURN
      IF(ITER.EQ.ITMAX) PAUSE 'Powell exceeding maximum iterations.'
      DO 14 J=1,N
        PTT(J)=2.d0*P(J)-PT(J)
        XIT(J)=P(J)-PT(J)
        PT(J)=P(J)
   14 CONTINUE
      FPTT=FUNC(PTT)
      IF(FPTT.GE.FP)GO TO 1
      T=2.d0*(FP-2.d0*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
      IF(T.GE.0.d0)GO TO 1
      CALL LINMIN(P,XIT,N,FRET)
      DO 15 J=1,N
        XI(J,IBIG)=XIT(J)
   15 CONTINUE
      GO TO 1
      END
      SUBROUTINE LINMIN(P,XI,N,FRET)
C Passed Parameters
      implicit none
      integer n
      double precision p(1),xi(1),fret
C Local Parameters
      integer j,ncom,nmax
      parameter (nmax=50)
      double precision ax,brent,bx,fa,fb,f1dim,func,fx,
     .  pcom,tol,xicom,xmin,xx
      parameter (tol=1.d-4)
      external f1dim
      common /f1com/ pcom(nmax),xicom(nmax),ncom
      NCOM=N
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
   11 CONTINUE
      AX=0
      XX=1
      BX=2
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
   12 CONTINUE
      RETURN
      END
      double precision FUNCTION F1DIM(X)
      implicit none
      double precision x
      integer nmax,ncom,j
      parameter (nmax=50)
      double precision pcom,xicom,xt(nmax),func
      common /f1com/ pcom(nmax),xicom(nmax),ncom

      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
   11 CONTINUE
      F1DIM=FUNC(XT)
C     print *, 'xt=', (xt(j), j=1,ncom), 'f=',f1dim
      RETURN
      END
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
      implicit none
C Passed Parameters
      double precision ax,bx,cx,fa,fb,fc,func
C Apollo bug requires that this be commented out:
      external func
C Local Parameters
      double precision fu,glimit,gold,q,r,tiny,u,ulim,dum
      PARAMETER (GLIMIT=100d0, TINY=1d-20)
      gold = (dsqrt(5d0)+1)/2
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
    1 IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2d0*dSIGN(MAX(dABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0d0)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0d0)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0d0)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      RETURN
      END
      double precision FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
C Passed Parameters
      implicit none
      double precision ax,bx,cx,f,tol,xmin
      external f
C Local Parameters
      double precision a,b,cgold,d,e,etemp,fu,fv,fw,fx
      double precision p,q,r,tol1,tol2,u,v,w,x,xm,zeps
      integer iter,itmax
      parameter (itmax=100,zeps=1d-10)
      cgold = 2/(3+dsqrt(5d0))
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0d0
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5d0*(A+B)
        TOL1=TOL*dABS(X)+ZEPS
        TOL2=2d0*TOL1
        IF(dABS(X-XM).LE.(TOL2-.5d0*(B-A))) GOTO 3
        IF(dABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2d0*(Q-R)
          IF(Q.GT.0d0) P=-P
          Q=dABS(Q)
          ETEMP=E
          E=D
          IF(dABS(P).GE.ABS(.5d0*Q*ETEMP).OR.P.LE.Q*(A-X).OR.
     *      P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=dSIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
    1   IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
    2   IF(dABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+dSIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
   11 CONTINUE
      PAUSE 'Brent exceed maximum iterations.'
    3 XMIN=X
      BRENT=FX
      RETURN
      END
