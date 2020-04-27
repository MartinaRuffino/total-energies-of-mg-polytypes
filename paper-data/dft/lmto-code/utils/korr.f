      SUBROUTINE KORR
     1(IELTER,BKOMMA,NACHKO,IREKOM,BKORRL,KONVKR,IFALLK,
     2TGRENZ,EPSILO,DELTAS,DELTAI,DELTAP,N,M,NS,NP,NY,
     3ZSTERN,XSTERN,ZBEST,X,S,P,Y,ZIELFU,RESTRI,GAUSSN,
     4GLEICH,TKONTR,KANAL)
C- H P Schwefel's evolution algorithm
C ----------------------------------------------------------------------
Ci Inputs:
Ci   IELTER: number of parents in a generation
Ci   BKOMMA: if F selection criterion applied to parents and offspring
Ci           if T selection criterion applied only to offspring
Ci   NACHKO: number of offspring from each generation
Ci           if BKOMMA=F, >=1; if T, > 6*IELTER
Ci   IREKOM: type of recombination
Ci           1, no recombination
Ci           2, discrete recombination of pairs of parents
Ci           3, Intermediate recombination of pairs of parents
Ci           4, Discrete recombination of all parents
Ci           5, Intermediate recombination of all parents in pairs
Ci   BKORRL: variability of mutation ellipsoid,
Ci           F, the ellipsoid cannot rotate
Ci           T, the ellipsoid can extend and rotate
Ci   KONVKR: type of convergence,
Ci           1, The difference in the objective function values between
Ci              the best and worst worst parents at the start of each
Ci              generation is used to determine whether to terminate the
Ci              search before the time limit is reached. Needs IELTER>1
Ci           2, (best >=2N) The change in the mean of all parental
Ci              objective function values in KONVKR generations is used
Ci              as the search termination criterion
Ci   TGRENZ: maximum time in seconds before run time exit
Ci   EPSILO: parameters to determine accuracy of optimisation at exit
Ci         (1) Lower bound to step sizes, absolute
Ci         (2) Lower bound to step sizes relative to values of
Ci             variables (not implemented)
Ci         (3) Limit to absolute value of objective function differences
Ci         (4) -ditto- but relative, not absolute
Ci   DELTAS: Factor used in step size change, hard wire to 1
Ci   DELTAI: similar to DELTAI, hard wire to 1
Ci   DELTAP: std deviation in random variation of the position angles
Ci           for the mutation ellipsoid, hard wire to 0.087
Ci        N: number of parameters
Ci        M: number of constraints
Ci       NS: length of array S, 1<= NS <= N
Ci       NP: length of array P, N*(NS-1)-((NS-1)*NS)/2 [BKORRL=T]
Ci       NY: length of array Y, (N+NS+NP+1)*IELTER*2   [BKORRL=T]
Ci   XSTERN: on input starting values of N parameters
Ci        X: work array of length N
Ci        S: NS step sizes, if too small may trap a local minimum
Ci        P: NP positional angles of ellipsoid, faut de mieux use 0.5
Ci        Y: work array length NY
Ci   ZIELFU: user supplied objective function: here f(x,n,ir) see below
Ci   RESTRI: user supplied constraints function. See remarks
Ci   GAUSSN: actually not passed as GAUSSN is built into the code below
Ci   GLEICH: user supplied function returning a double precision
Ci           random number in [0,1]
Ci   TKONTR: user supplied function to return elapsed CPU time in seconds
Ci    KANAL: logical unit for standard out
Co Outputs:
Co   IFALLK: return code, normal exit is 2
Co   ZSTERN: best value of (constrained) objective function
Ci   XSTERN: on output optimised values of N parameters
Co    ZBEST: current best value of objective function
Cr Remarks
Cr  This is the FORTRAN77 code for H P Schwefel's evolution algorithm
Cr  practically unchanged except for the calls to the objective
Cr  function (Zielfunktion) to align with the praxis-like calls in fmin.
Cr
Cr  The constraints function has the form
Cr        double precision function restri(j,n,x)
Cr        integer j,n
Cr        double precision x(n)
Cr  and returns restri < 0 if the j'th constraint is not satisfied
Cr
Cr Reference
Cr  "Numerical Optimization of Computer Models," H-P Schwefel, John
Cr  Wiley, 1981, ISBN-10: 0471099880, ISBN-13: 978-0471099888
C ----------------------------------------------------------------------
      implicit none
      LOGICAL BKOMMA,BKORRL,BFATAL,BKONVG,BLETAL
      double precision EPSILO(4),XSTERN(N),X(N),S(NS),P(NP),
     1Y(NY)
      double precision PIHALB,PIEINS,PI3HLB,PIZWEI,D
      double precision ZSTERN,TGRENZ,DELTAS,DELTAI,DELTAP
      double precision ZBEST,TMAXIM,Z1,DSMAXI,DPMAXI,Z,ZSCHL,Z2
      double precision CM1,C1
      integer IELTER,NACHKO,IREKOM,KONVKR,IFALLK
      integer Izero,N,M,NS,NP,NY,KANAL,NL,NM,NZ,LBEST
      integer L,I,L1,L2,KONVZ,L3,LMUTAT,LSCHL,K1,K2,K
      integer ir
      COMMON/PIDATA/PIHALB,PIEINS,PI3HLB,PIZWEI
      double precision,external :: RESTRI,GAUSSN,GLEICH
      double precision,external :: ZULASS,TKONTR
      double precision ZIELFU
      integer iprint,i1mach
      Izero = 0
      ir = 0
      CM1   = -1.0d0
      C1    =  1.0d0
      CALL PRUEFG
     1(IELTER,BKOMMA,NACHKO,IREKOM,BKORRL,KONVKR,TGRENZ,
     2EPSILO,DELTAS,DELTAI,DELTAP,N,M,NS,NP,NY,KANAL,
     3BFATAL)
C
C   CHECK INPUT PARAMETERS FOR FORMAL ERRORS.
C
      IF(BFATAL) RETURN
C
C   PREPARE AUXILIARY  QUANTITIES. TIMING  MONITORED IN
C   ACCORDANCE  WITH  THE  TKONTR  FUNCTION  FROM  HERE
C   ONWARDS.
C
      TMAXIM=TGRENZ+TKONTR(D)
      IF(.NOT.BKORRL) GOTO 1
      PIHALB=2.*ATAN(1.)
      PIEINS=PIHALB+PIHALB
      PI3HLB=PIEINS+PIHALB
      PIZWEI=PIEINS+PIEINS
      NL=1+N-NS
      NM=N-1
  1   NZ=NY/(IELTER+IELTER)
      IF(M.EQ.0) GOTO 2
C
C   CHECK FEASIBILITY OF INITIAL VECTOR XSTERN.
C
      IFALLK=-1
      ZSTERN=ZULASS(N,M,XSTERN,RESTRI)
      IF(ZSTERN.GT.0.) GOTO 3
  2   IFALLK=1
C      ZSTERN=ZIELFU(N,XSTERN)
      ZSTERN=ZIELFU(XSTERN,N,ir)

      if (iprint() >= 30) then
        call awrit3('KORR: Evaluation of objective function for '//
     .    'parent 1.%N objective function, Z*=%d%N'//
     .    ' parameters, X*=%n:1d',' ',2056,i1mach(2),ZSTERN,N,XSTERN)
      endif

  3   CALL SPEICH
     1(Izero,BKORRL,EPSILO,N,NS,NP,NY,ZSTERN,XSTERN,S,P,Y)
C
C   THE INITIAL VALUES SUPPLIED  BY THE USER ARE STORED
C   IN FIELD Y AS THE DATA OF THE FIRST PARENT.
C
      IF(KONVKR.GT.1) Z1=ZSTERN
      ZBEST=ZSTERN
      LBEST=0
      IF(IELTER.EQ.1) GOTO 16
      DSMAXI=DELTAS
      DPMAXI=MIN(DELTAP*10.0d0,PIHALB)
      DO 14 L=2,IELTER
C
C   IF IELTER > 1, THE OTHER IELTER - 1 INITIAL VECTORS
C   ARE DERIVED FROM THE VECTOR FOR THE FIRST PARENT BY
C   MUTATION (WITHOUT SELECTION). THE STRATEGY
C   PARAMETERS ARE WIDELY SPREAD.
C
      DO 4 I=1,NS
  4   S(I)=Y(N+I)
  5   IF(.NOT.BKORRL) GOTO 7
      DO 6 I=1,NP
  6   P(I)=Y(N+NS+I)
  7   CALL MUTATI
     1(NL,NM,BKORRL,DSMAXI,DELTAI,DPMAXI,N,NS,NP,X,S,P,
     2GAUSSN,GLEICH)
C
C   MUTATION IN ALL OBJECT AND STRATEGY PARAMETERS.
C
      DO 8 I=1,N
  8   X(I)=X(I)+Y(I)
      IF(IFALLK.GT.0) GOTO 9
C
C   IF THE STARTING POINT IS  INFEASIBLE, EACH MUTATION
C   IS CHECKED AT ONCE TO SEE WHETHER A FEASIBLE VECTOR
C   HAS BEEN FOUND.  THE SEARCH ENDS WITH IFALLK = 0 IF
C   THIS IS SO.
C
      Z=ZULASS(N,M,X,RESTRI)
      IF(Z)40,40,12
  9   IF(M.EQ.0) GOTO 11
      IF(.NOT.BLETAL(N,M,X,RESTRI)) GOTO 11
C
C   IF  A  MUTATION  FROM  A  FEASIBLE  STARTING  POINT
C   RESULTS IN  AN INFEASIBLE  X VECTOR, THEN  THE STEP
C   SIZES ARE REDUCED (ON THE ASSUMPTION THAT THEY WERE
C   INITIALLY  TOO   LARGE)  IN  ORDER   TO  AVOID  THE
C   THE COMSUMPTION  OF EXCESSIVE TIME  IN DEFINING THE
C   THE FIRST PARENT GENERATION.
C
      DO 10 I=1,NS
 10   S(I)=S(I)*.5
      GOTO 5
C 11   Z=ZIELFU(N,X)
 11   Z=ZIELFU(X,N,ir)

      if (iprint() >= 30) then
        call awrit4('KORR: Evaluation of objective function for '//
     .    'parent %i.%N objective function, Z=%d%N'//
     .    ' parameters, X=%n:1d',' ',2056,i1mach(2),L,Z,N,X)
      endif

 12   IF(Z.GT.ZBEST) GOTO 13
      ZBEST=Z
      LBEST=L-1
      DSMAXI=DSMAXI*LOG(2.0d0)
 13   CALL SPEICH
     1((L-1)*NZ,BKORRL,EPSILO,N,NS,NP,NY,Z,X,S,P,Y)
C
C   STORE PARENT DATA IN ARRAY Y.
C
      IF(KONVKR.GT.1) Z1=Z1+Z
 14   CONTINUE
C
C   THE  INITIAL  PARENT GENERATION  IS  NOW  COMPLETE.
C   ZSTERN AND XSTERN,  WHICH HOLD THE BEST VALUES, ARE
C   OVERWRITTEN  WHEN  AN IMPROVEMENT  OF  THE  INITIAL
C   SITUATION IS OBTAINED.
C
      IF(LBEST.EQ.0) GOTO 16
      ZSTERN=ZBEST
      K=LBEST*NZ
      DO 15 I=1,N
 15   XSTERN(I)=Y(K+I)
 16   L1=IELTER
      L2=0
      IF(KONVKR.GT.1) KONVZ=0

      if (iprint() >= 30) then
        call awrit3('KORR: Initial parent generation complete.%N'//
     .    ' objective function, Z*=%d%N parameters, X*=%n:1d',' ',2056,
     .    i1mach(2),ZSTERN,N,XSTERN)
      endif
C
C   ALL INITIALIZATION STEPS COMPLETED AT THIS POINT.
C   EACH FRESH GENERATION NOW STARTS AT LABEL 17.
C
 17   L3=L2
      L2=L1
      L1=L3
      IF(M.GT.0) L3=0
      LMUTAT=0
C
C   LMUTAT IS THE MUTATION COUNTER WITHIN A GENERATION,
C   WHILE L3 IS  THE COUNTER FOR LETHAL  MUTATIONS WHEN
C   THE PROBLEM INVOLVES CONSTRAINTS.
C
      IF(BKOMMA) GOTO 18
C
C   IF BKOMMA=.FALSE.  HAS BEEN  SELECTED,  THE PARENTS
C   MUST BE INCORPORATED IN THE SELECTION. THE DATA FOR
C   THESE ARE  TRANSFERRED FROM  THE FIRST  (OR SECOND)
C   PART OF THE ARRAY Y TO  THE SECOND (OR FIRST) PART.
C   IN THIS  CASE  THE  WORST INDIVIDUAL  MUST  ALSO BE
C   KNOWN,  THIS  IS   REPLACED  BY  THE  FIRST  BETTER
C   DESCENDANT.
C
      CALL UMSPEI
     1(L1*NZ,L2*NZ,IELTER*NZ,NY,Y)
      CALL MINMAX
     1(CM1,L2,NZ,ZSCHL,LSCHL,IELTER,NY,Y)
C
C   THE GENERATION OF EACH DESCENDANT STARTS AT LABEL 18
C
 18   IF(IREKOM.GT.3) GOTO 19
C
C   RANDOM CHOICE  OF A PARENT OR  OF A PAIR OF PARENTS
C   IN ACCORDANCE WITH THE  VALUE CHOSEN FOR IREKOM. IF
C   IREKOM=3 OR IREKOM=5, THE CHOICE OF PARENTS IS MADE
C   WITHIN GNPOOL.
C
      K1=L1+IELTER*GLEICH(D)
      IF(IREKOM.GT.1) K2=L1+IELTER*GLEICH(D)
 19   CALL GNPOOL
     1(1,L1,K1,K2,NZ,N,IELTER,IREKOM,NS,NY,S,Y,GLEICH)
C
C   STEP SIZES  SUPPLIED  FOR  THE DESCENDANT  FROM THE
C   POOL OF GENES.
C
      IF(BKORRL) CALL GNPOOL
     1(2,L1,K1,K2,NZ,N+NS,IELTER,IREKOM,NP,NY,P,Y,GLEICH)
C
C   POSITIONAL  ANGLES OF  ELLIPSOID  SUPPLIED  FOR THE
C   DESCENDANT FROM THE POOL  OF GENES WHEN CORRELATION
C   IS REQUIRED.
C
      CALL MUTATI
     1(NL,NM,BKORRL,DELTAS,DELTAI,DELTAP,N,NS,NP,X,S,P,
     2GAUSSN,GLEICH)
C
C   CALL TO  MUTATION  SUBROUTINE  FOR  ALL  VARIABLES,
C   INCLUDING  POSSIBLY  COORDINATE  TRANSFORMATION.  S
C   (AND P)  ARE  ALREADY  THE  NEW  ATTRIBUTES  OF THE
C   DESCENDANT,  WHILE X REPRESENTS  THE CHANGES  TO BE
C   MADE IN THE OBJECT VARIABLES.
C
      CALL GNPOOL
     1(3,L1,K1,K2,NZ,0,IELTER,IREKOM,N,NY,X,Y,GLEICH)
C
C   OBJECT VARIABLES  SUPPLIED FOR  THE DESCENDANT FROM
C   THE POOL OF GENES  AND ADDITION OF THE MODIFICATION
C   VECTOR. X NOW  REPRESENTS  THE  NEW  STATE  OF  THE
C   DESCENDANT.
C
      LMUTAT=LMUTAT+1
      IF(IFALLK.GT.0) GOTO 20
C
C   EVALUATION OF THE AUXILIARY  OBJECTIVE FUNCTION FOR
C   THE SEARCH FOR A FEASIBLE VECTOR.
C
      Z=ZULASS(N,M,X,RESTRI)
      IF(Z)40,40,22
 20   IF(M.EQ.0) GOTO 21
C
C   CHECK FEASIBILITY  OF DESCENDANT.  IF THE RESULT IS
C   NEGATIVE  (LETHAL MUTATION),  THE  MUTATION  IS NOT
C   COUNTED AS REGARDS THE NACHKO PARAMETER.
C
      IF(.NOT.BLETAL(N,M,X,RESTRI)) GOTO 21
      IF(.NOT.BKOMMA) GOTO 25
      LMUTAT=LMUTAT-1
      L3=L3+1
      IF(L3.LT.NACHKO) GOTO 18
      L3=0
C
C   TIME CHECK MADE NOT ONLY  AFTER EACH GENERATION BUT
C   ALSO  AFTER  EVERY  NACHKO   LETHAL  MUTATIONS  FOR
C   CERTAINTY.
C
      IF(TKONTR(D).LT.TMAXIM) GOTO 18
      IFALLK=3
      GOTO 26
C 21   Z=ZIELFU(N,X)
 21   Z=ZIELFU(X,N,ir)
C
C   EVALUATION  OF  OBJECTIVE  FUNCTION  VALUE  FOR THE
C   DESCENDANT.
C
      if (iprint() >= 30) then
        call awrit4('KORR: Evaluation of objective function for '//
     .    'descendent %i.%N objective function, Z=%d%N'//
     .    ' parameters, X=%n:1d',' ',2056,i1mach(2),LMUTAT,Z,N,X)
      endif

 22   IF(BKOMMA.AND.LMUTAT.LE.IELTER) GOTO 23
      IF(Z-ZSCHL)24,24,25
 23   LSCHL=L2+LMUTAT-1
 24   CALL SPEICH
     1(LSCHL*NZ,BKORRL,EPSILO,N,NS,NP,NY,Z,X,S,P,Y)
C
C   TRANSFER OF DATA  OF DESCENDANT  TO PART OF ARRAY Y
C   HOLDING THE PARENTS FOR THE NEXT GENERATION.
C
      IF(.NOT.BKOMMA.OR.LMUTAT.GE.IELTER) CALL MINMAX
     1(CM1,L2,NZ,ZSCHL,LSCHL,IELTER,NY,Y)
C
C   LOOK FOR THE  CURRENTLY WORST  INDIVIDUAL STORED IN
C   ARRAY Y WITHOUT CONSIDERING  THE PARENTS THAT STILL
C   CAN PRODUCE DESCENDANTS IN THIS GENERATION.
C
 25   IF(LMUTAT.LT.NACHKO) GOTO 18
C
C   END OF GENERATION.
C
 26   CALL MINMAX
     1(C1,L2,NZ,ZBEST,LBEST,IELTER,NY,Y)
C
C   LOOK FOR  THE  BEST  OF  THE  INDIVIDUALS  HELD  AS
C   PARENTS FOR THE NEXT  GENERATION. IF THIS IS BETTER
C   THAN ANY DESCENDANT  PREVIOUSLY GENERATED, THE DATA
C   ARE WRITTEN INTO ZSTERN AND XSTERN.
C
      IF(ZBEST.GT.ZSTERN) GOTO 28
      ZSTERN=ZBEST
      K=LBEST*NZ
      DO 27 I=1,N
 27   XSTERN(I)=Y(K+I)
 28   continue
      if (iprint() >= 30) then
        call awrit3('KORR: generation complete.%N'//
     .    ' objective function, Z*=%d%N parameters, X*=%n:1d',' ',2056,
     .    i1mach(2),ZSTERN,N,XSTERN)
      endif

      IF(IFALLK.EQ.3) GOTO 30
      Z2=0.0d0
      K=L2*NZ
      DO 29 L=1,IELTER
      K=K+NZ
 29   Z2=Z2+Y(K)
      CALL ABSCHA
     1(IELTER,KONVKR,IFALLK,EPSILO,ZBEST,ZSCHL,Z1,Z2,
     2KONVZ,BKONVG)
C
C   TEST CONVERGENCE CRITERION.
C
      IF(BKONVG) GOTO 30
C
C   CHECK TIME ELAPSED.
C
      IF(TKONTR(D).LT.TMAXIM) GOTO 17
C
C   PREPARE  FINAL DATA  FOR  RETURN  FROM KORR  IF THE
C   STARTING POINT WAS FEASIBLE.
C
 30   K=LBEST*NZ
      DO 31 I=1,N
      K=K+1
 31   X(I)=Y(K)
      DO 32 I=1,NS
      K=K+1
 32   S(I)=Y(K)
      IF(.NOT.BKORRL) RETURN
      DO 33 I=1,NP
      K=K+1
 33   P(I)=Y(K)
      RETURN
C
C   PREPARE  FINAL DATA  FOR  RETURN  FROM KORR  IF THE
C   STARTING POINT WAS INFEASIBLE.
C
 40   DO 41 I=1,N
 41   XSTERN(I)=X(I)
C      ZSTERN=ZIELFU(N,XSTERN)
      ZSTERN=ZIELFU(XSTERN,N,ir)
      ZBEST=ZSTERN
      IFALLK=0
      RETURN
      END
      SUBROUTINE ABSCHA
     1(IELTER,KONVKR,IFALLK,EPSILO,ZBEST,ZSCHL,Z1,Z2,
     2KONVZ,BKONVG)
      implicit none
      LOGICAL BKONVG
      double precision EPSILO(4)
      double precision Z1,Z2,ZSCHL,ZBEST,DELTAF
      integer KONVKR,KONVZ,IELTER,IFALLK
      integer iprint,i1mach
      IF(KONVKR.EQ.1) GOTO 1
      KONVZ=KONVZ+1

      if (iprint() >= 30) then
      call awrit2(' ABSCHA: convergence test, KONVKR=%i, KONVZ=%i',' ',
     .    128,i1mach(2),KONVKR,KONVZ)
      endif

      IF(KONVZ.LT.KONVKR) GOTO 3
      KONVZ=0
      DELTAF=Z1-Z2
      Z1=Z2

      if (iprint() >= 30) then
      call awrit6(' ABSCHA: convergence test, KONVKR=%i, KONVZ=%i.'/
     .    /' Z2=%d'/
     .    /' DELTAF=%d. Compare with eps(3)*N=%d, eps(4)*|Z2|=%d.',' ',
     .    1028,i1mach(2),KONVKR,KONVZ,Z2,DELTAF,EPSILO(3)*IELTER,
     .    EPSILO(4)*ABS(Z2))
      endif

      GOTO 2
1     DELTAF=(ZSCHL-ZBEST)*IELTER

      if (iprint() >= 30) then
      call awrit6(' ABSCHA: convergence test, KONVKR=%i. Z1=%d, Z2=%d,'/
     .    /' DELTAF=%d. Compare with eps(3)*N=%d, eps(4)*|Z2|=%d.',' ',
     .    1028,i1mach(2),KONVKR,Z1,Z2,DELTAF,EPSILO(3)*IELTER,
     .    EPSILO(4)*ABS(Z2))
      endif

2     IF(DELTAF.GT.EPSILO(3)*IELTER) GOTO 3
      IF(DELTAF.GT.EPSILO(4)*ABS(Z2)) GOTO 3
      IFALLK=ISIGN(2,IFALLK)
      BKONVG=.TRUE.
      RETURN
3     BKONVG=.FALSE.
      RETURN
      END
      LOGICAL FUNCTION BLETAL
     1(N,M,X,RESTRI)
      implicit none
      double precision X(N)
      integer J,M,N
      double precision RESTRI
      external RESTRI
      DO 1 J=1,M
      IF(RESTRI(J,N,X).LT.0.) GOTO 2
1     CONTINUE
      BLETAL=.FALSE.
      RETURN
2     BLETAL=.TRUE.
      RETURN
      END
      SUBROUTINE DREHNG
     1(NL,NM,N,NP,X,P)
      implicit none
      double precision X(N),P(NP)
      double precision CO,X1,X2,SI
      integer N,NP,NQ,NL,NM,N1,N2,I,II
      NQ=NP
      DO 1 II=NL,NM
      N1=N-II
      N2=N
      DO 1 I=1,II
      X1=X(N1)
      X2=X(N2)
      SI=SIN(P(NQ))
      CO=COS(P(NQ))
      X(N2)=X1*SI+X2*CO
      X(N1)=X1*CO-X2*SI
      N2=N2-1
1     NQ=NQ-1
      RETURN
      END
      SUBROUTINE EVOL(N,M,LF,LR,LS,TM,EA,EB,EC,ED,SN,FB,
     1XB,SM,X,F,G,T,Z,R)
      DIMENSION XB(1),SM(1),X(1),L(10)
      COMMON/EVZ/LZ
      EXTERNAL R
      TN=TM+T(D)
      LZ=1
      IF(M)4,4,1
1     LF=-1
      FB=0.
      DO 3 J=1,M
      FG=G(J,N,XB)
      IF(FG)2,3,3
2     FB=FB-FG
3     CONTINUE
      IF(FB)4,4,5
4     LF=1
      FB=F(N,XB)
5     DO 6 K=1,10
6     L(K)=N*K/5
      LE=N+N
      LM=0
      LC=0
      FC=FB
7     DO 8 I=1,N
8     X(I)=XB(I)+Z(SM(I),R)
      IF(LF)9,9,12
9     FF=0.
      DO 11 J=1,M
      FG=G(J,N,X)
      IF(FG)10,11,11
10    FF=FF-FG
11    CONTINUE
      IF(FF)32,32,16
12    IF(M)15,15,13
13    DO 14 J=1,M
      IF(G(J,N,X))19,14,14
14    CONTINUE
15    FF=F(N,X)
16    IF(FF-FB)17,17,19
17    LE=LE+1
      FB=FF
      DO 18 I=1,N
18    XB(I)=X(I)
19    LM=LM+1
      IF(LM-N*LR)7,20,20
20    K=1
      IF(LE-L(1)-N-N)23,22,21
21    K=K-1
22    K=K-1
23    DO 24 I=1,N
24    SM(I)=AMAX1(SM(I)*SN**K,ABS(XB(I))*EB,EA)
      DO 25 K=1,9
25    L(K)=L(K+1)
      L(10)=LE
      LM=0
      LC=LC+1
      IF(LC-10*LS)31,26,26
26    IF(FC-FB-EC)28,28,27
27    IF((FC-FB)/ED-ABS(FC))28,28,30
28    LF=ISIGN(2,LF)
29    RETURN
30    LC=0
      FC=FB
31    IF(T(D)-TN)7,29,29
32    DO 33 I=1,N
33    XB(I)=X(I)
      FB=F(N,XB)
      LF=0
      GOTO 29
      END
      DOUBLE PRECISION FUNCTION GAUSSN
     1(SIGMA,GLEICH)
      implicit none
      double precision GLEICH
      double precision SIGMA,D,U,U0,U1,U2,X,Y
1     U=GLEICH(D)
      U0=GLEICH(D)
      IF(U.GE..919544406) GOTO 2
      X=2.40375766*(U0+U*.825339283)-2.11402808
      GOTO 10
2     IF(U.LT..965487131) GOTO 4
3     U1=GLEICH(D)
      Y=SQRT(4.46911474-2.*LOG(U1))
      U2=GLEICH(D)
      IF(Y*U2.GT.2.11402808) GOTO 3
      GOTO 9
4     IF(U.LT..949990709) GOTO 6
5     U1=GLEICH(D)
      Y=1.84039875+U1*.273629336
      U2=GLEICH(D)
      IF(.398942280*EXP(-.5*Y*Y)-.443299126+Y*.209694057
     1.LT.U2*.0427025816) GOTO 5
      GOTO 9
6     IF(U.LT..925852334) GOTO 8
7     U1=GLEICH(D)
      Y=.289729574+U1*1.55066917
      U2=GLEICH(D)
      IF(.398942280*EXP(-.5*Y*Y)-.443299126+Y*.209694057
     1.LT.U2*.0159745227) GOTO 7
      GOTO 9
8     U1=GLEICH(D)
      Y=U1*.289729574
      U2=GLEICH(D)
      IF(.398942280*EXP(-.5*Y*Y)-.382544556
     1.LT.U2*.0163977244) GOTO 8
9     X=Y
      IF(U0.GE..5) X=-Y
10    GAUSSN=SIGMA*X
      RETURN
      END
      SUBROUTINE GNPOOL
     1(J,L1,K1,K2,NZ,NN,IELTER,IREKOM,NX,NY,XX,Y,GLEICH)
      implicit none
      double precision  XX(NX),Y(NY)
      double precision PIHALB,PIEINS,PI3HLB,PIZWEI
      double precision D,XX1,XX2,XXI,DXX,ADD
      integer J,L1,K1,K2,NZ,NN,IELTER,IREKOM,NX,NY
      integer KI1,I,KI2,KI
      COMMON/PIDATA/PIHALB,PIEINS,PI3HLB,PIZWEI
      double precision,external :: GLEICH
      IF(J.EQ.3) GOTO 11
      GOTO(1,1,1,7,9),IREKOM
1     KI1=K1*NZ+NN
      IF(IREKOM.GT.1) GOTO 3
      DO 2 I=1,NX
2     XX(I)=Y(KI1+I)
      RETURN
3     KI2=K2*NZ+NN
      IF(IREKOM.EQ.3) GOTO 5
      DO 4 I=1,NX
      KI=KI1
      IF(GLEICH(D).GE..5) KI=KI2
4     XX(I)=Y(KI+I)
      RETURN
5     DO 6 I=1,NX
      XX1=Y(KI1+I)
      XX2=Y(KI2+I)
      XXI=(XX1+XX2)*.5
      IF(J.EQ.1) GOTO 6
      DXX=XX1-XX2
      IF(ABS(DXX).LT.PIHALB) GOTO 6
      ADD=SIGN(PIHALB,DXX)
      XXI=XXI+ADD
      IF(ABS(DXX).GE.PI3HLB) XXI=XXI+ADD
6     XX(I)=XXI
      RETURN
7     DO 8 I=1,NX
8     XX(I)=Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I)
      RETURN
9     DO 10 I=1,NX
      XX1=Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I)
      XX2=Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I)
      XXI=(XX1+XX2)*.5
      IF(J.EQ.1) GOTO 10
      DXX=XX1-XX2
      IF(ABS(DXX).LT.PIHALB) GOTO 10
      ADD=SIGN(PIHALB,DXX)
      XXI=XXI+ADD
      IF(ABS(DXX).GE.PI3HLB) XXI=XXI+ADD
10    XX(I)=XXI
      RETURN
11    GOTO(12,12,12,18,20),IREKOM
12    KI1=K1*NZ+NN
      IF(IREKOM.GT.1) GOTO 14
      DO 13 I=1,NX
13    XX(I)=XX(I)+Y(KI1+I)
      RETURN
14    KI2=K2*NZ+NN
      IF(IREKOM.EQ.3) GOTO 16
      DO 15 I=1,NX
      KI=KI1
      IF(GLEICH(D).GE..5) KI=KI2
15    XX(I)=XX(I)+Y(KI+I)
      RETURN
16    DO 17 I=1,NX
17    XX(I)=XX(I)+(Y(KI1+I)+Y(KI2+I))*.5
      RETURN
18    DO 19 I=1,NX
19    XX(I)=XX(I)+Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I)
      RETURN
20    DO 21 I=1,NX
21    XX(I)=XX(I)+(Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I)
     1+Y((L1+INT(IELTER*GLEICH(D)))*NZ+NN+I))*.5
      RETURN
      END
      SUBROUTINE MINMAX
     1(C,LL,NZ,ZM,LM,IELTER,NY,Y)
      implicit none
      double precision C,ZM,ZZ
      double precision Y(NY)
      integer LL,NZ,LM,IELTER,NY,K1,K2,KM,K
      LM=LL
      K1=LL*NZ+NZ
      ZM=Y(K1)
      IF(IELTER.EQ.1) RETURN
      K1=K1+NZ
      K2=(LL+IELTER)*NZ
      KM=LL
      DO 1 K=K1,K2,NZ
      KM=KM+1
      ZZ=Y(K)
      IF((ZZ-ZM)*C.GT.0.) GOTO 1
      ZM=ZZ
      LM=KM
1     CONTINUE
      RETURN
      END
      SUBROUTINE MUTATI
     1(NL,NM,BKORRL,DELTAS,DELTAI,DELTAP,N,NS,NP,X,S,P,
     2GAUSSN,GLEICH)
      implicit none
      LOGICAL BKORRL
      double precision DELTAS,DELTAI,DELTAP
      double precision X(N),S(NS),P(NP),DS
      double precision GAUSSN
      integer NL,NM,N,NS,NP,I
      double precision,external :: GLEICH
      DS=GAUSSN(DELTAS,GLEICH)
      DO 1 I=1,NS
1     S(I)=S(I)*EXP(DS+GAUSSN(DELTAI,GLEICH))
      DO 2 I=1,N
2     X(I)=GAUSSN(S(MIN0(I,NS)),GLEICH)
      IF(.NOT.BKORRL) RETURN
      DO 3 I=1,NP
3     P(I)=P(I)+GAUSSN(DELTAP,GLEICH)
      CALL DREHNG
     1(NL,NM,N,NP,X,P)
      RETURN
      END
      SUBROUTINE PRUEFG
     1(IELTER,BKOMMA,NACHKO,IREKOM,BKORRL,KONVKR,TGRENZ,
     2EPSILO,DELTAS,DELTAI,DELTAP,N,M,NS,NP,NY,KANAL,
     3BFATAL)
      LOGICAL BKOMMA,BKORRL,BFATAL
      double precision EPSILO(4),TGRENZ,DELTAS,DELTAI,DELTAP
      integer IELTER, NACHKO,IREKOM,KONVKR,N,M,NS,NP,NY
      integer KANAL
100   FORMAT(1H ,' CORRECTION. IELTER >0 . ASSUMED: 2 AND
     1 KONVKR = ',I5)
101   FORMAT(1H ,' CORRECTION. NACHKO >0 . ASSUMED: ',I5)
102   FORMAT(1H ,'WARNING.BETTER VALUE NACHKO>=6*IELTER')
103   FORMAT(1H ,' CORRECTION. IF BKOMMA = .TRUE., THEN
     1 NACHKO > IELTER . ASSUMED: ',I3)
104   FORMAT(1H ,' CORRECTION. 0< IREKOM <6 .ASSUMED: 1')
105   FORMAT(1H ,' CORRECTION. IF IELTER = 1, THEN
     1 IREKOM = 1 . ASSUMED: 1')
106   FORMAT(1H ,' CORRECTION. IF N = 1 OR NS = 1, THEN
     1 BKORRL = .FALSE. . ASSUMED: .FALSE.')
107   FORMAT(1H ,' CORRECTION. KONVKR >0 . ASSUMED: ',I5)
108   FORMAT(1H ,' CORRECTION. IF IELTER = 1, THEN
     1 KONVKR > 1 . ASSUMED: ',I5)
109   FORMAT(1H ,' CORRECTION. EPSILO(',I1,') > 0. .
     1 SIGN REVERSED')
110   FORMAT(1H ,' WARNING. EPSILO(',I1,') TOO SMALL.
     1 TREATED AS 0. .')
111   FORMAT(1H ,' CORRECTION. DELTAS >= 0. .
     1 SIGN REVERSED')
112   FORMAT(1H ,' WARNING. EXP(DELTAS) = 1.
     1 OVER-ALL STEP SIZE CONSTANT')
113   FORMAT(1H ,' CORRECTION. DELTAI >= 0. .
     1 SIGN REVERSED')
114   FORMAT(1H ,' WARNING. EXP(DELTAI) = 1.
     1 STEP-SIZE RELATIONS CONSTANT')
115   FORMAT(1H ,' CORRECTION. DELTAP >= 0. .
     1 SIGN REVERSED')
116   FORMAT(1H ,' WARNING. DELTAP = 0.
     1 CORRELATION REMAINS FIXED')
117   FORMAT(1H ,' WARNING. TGRENZ <= 0.
     1 ONE GENERATION TESTED')
118   FORMAT(1H ,' CORRECTION. M >= 0 . ASSUMED: 0')
119   FORMAT(1H ,' FATAL ERROR. N <= 0')
120   FORMAT(1H ,' FATAL ERROR. NS <= 0')
121   FORMAT(1H ,' FATAL ERROR. NP <= 0')
122   FORMAT(1H ,' CORRECTION. 1<= NS <=N .ASSUMED: ',I5)
123   FORMAT(1H ,' CORRECTION. IF BKOORL = .FALSE., THEN
     1 NP = 1 . ASSUMED: 1')
124   FORMAT(1H ,' FATAL ERROR. NY < (N+NS+1)*IELTER*2')
125   FORMAT(1H ,' CORRECTION. NY = (N+NS+1)*IELTER*2 .
     1 ASSUMED: ',I5)
126   FORMAT(1H ,'FATAL ERROR.NP<N*(NS-1)-((NS-1)*NS)/2')
127   FORMAT(1H ,' CORRECTION. NP=N*(NS-1)-((NS-1)*NS)/2.
     1 ASSUMED: ',I5)
128   FORMAT(1H ,' FATAL ERROR. NY<(N+NS+NP+1)*IELTER*2')
129   FORMAT(1H ,' CORRECTION. NY = (N+NS+NP+1)*IELTER*2.
     1 ASSUMED: ',I5)
      BFATAL=.TRUE.
      IF(IELTER.GT.0) GOTO 1
      IELTER=2
      KONVKR=N+N
      WRITE(KANAL,100)KONVKR
1     IF(NACHKO.GT.0) GOTO 2
      NACHKO=6*IELTER
      WRITE(KANAL,101)NACHKO
2     IF(.NOT.BKOMMA.OR.NACHKO.GE.6*IELTER) GOTO 3
      WRITE(KANAL,102)
      IF(NACHKO.GT.IELTER) GOTO 3
      NACHKO=6*IELTER
      WRITE(KANAL,103)NACHKO
3     IF(IREKOM.GT.0.AND.IREKOM.LT.6) GOTO 4
      IREKOM=1
      WRITE(KANAL,104)
      GOTO 5
4     IF(IREKOM.EQ.1.OR.IELTER.NE.1) GOTO 5
      IREKOM=1
      WRITE(KANAL,105)
5     IF(.NOT.BKORRL.OR.(N.GT.1.AND.NS.GT.1)) GOTO 6
      BKORRL=.FALSE.
      WRITE(KANAL,106)
6     IF(KONVKR.GT.0) GOTO 7
      IF(IELTER.EQ.1) KONVKR=N+N
      IF(IELTER.GT.1) KONVKR=1
      WRITE(KANAL,107)KONVKR
      GOTO 8
7     IF(KONVKR.GT.1.OR.IELTER.GT.1) GOTO 8
      KONVKR=N+N
      WRITE(KANAL,108)KONVKR
8     DO 12 I=1,4
      IF(I.EQ.2.OR.I.EQ.4) GOTO 9
      IF(EPSILO(I))10,11,12
9     IF((1.+EPSILO(I))-1.)10,11,12
10    EPSILO(I)=-EPSILO(I)
      WRITE(KANAL,109)I
      GOTO 12
11    WRITE(KANAL,110)I
12    CONTINUE
      IF(EXP(DELTAS)-1.)13,14,15
13    DELTAS=-DELTAS
      WRITE(KANAL,111)
      GOTO 15
14    IF(EXP(DELTAI).NE.1.) GOTO 15
      WRITE(KANAL,112)
15    IF(EXP(DELTAI)-1.)16,17,18
16    DELTAI=-DELTAI
      WRITE(KANAL,113)
      GOTO 18
17    IF(IREKOM.GT.1.AND.EXP(DELTAS).GT.1.) GOTO 18
      WRITE(KANAL,114)
18    IF(.NOT.BKORRL) GOTO 21
      IF(DELTAP)19,20,21
19    DELTAP=-DELTAP
      WRITE(KANAL,115)
      GOTO 21
20    WRITE(KANAL,116)
21    IF(TGRENZ.GT.0.) GOTO 22
      WRITE(KANAL,117)
22    IF(M.GE.0) GOTO 23
      M=0
      WRITE(KANAL,118)
23    IF(N.GT.0) GOTO 24
      WRITE(KANAL,119)
      RETURN
24    IF(NS.GT.0) GOTO 25
      WRITE(KANAL,120)
      RETURN
25    IF(NP.GT.0) GOTO 26
      WRITE(KANAL,121)
      RETURN
26    IF(NS.LE.N) GOTO 27
      NS=N
      WRITE(KANAL,122)N
27    IF(BKORRL) GOTO 31
      IF(NP.EQ.1) GOTO 28
      NP=1
      WRITE(KANAL,123)
28    NYY=(N+NS+1)*IELTER*2
      IF(NY-NYY)29,37,30
29    WRITE(KANAL,124)
      RETURN
30    NY=NYY
      WRITE(KANAL,125)NY
      GOTO 37
31    NPP=N*(NS-1)-((NS-1)*NS)/2
      IF(NP-NPP)32,34,33
32    WRITE(KANAL,126)
      RETURN
33    NP=NPP
      WRITE(KANAL,127)NP
34    NYY=(N+NS+NP+1)*IELTER*2
      IF(NY-NYY)35,37,36
35    WRITE(KANAL,128)
      RETURN
36    NY=NYY
      WRITE(KANAL,129)NY
37    BFATAL=.FALSE.
      RETURN
      END
      SUBROUTINE SPEICH
     1(J,BKORRL,EPSILO,N,NS,NP,NY,ZZ,XX,S,P,Y)
      implicit none
      LOGICAL BKORRL
      double precision EPSILO(4),XX(N),S(NS),P(NP),Y(NY)
      double precision PIHALB,PIEINS,PI3HLB,PIZWEI,PI,ZZ
      integer K,J,I,N,NS,NP,NY
      COMMON/PIDATA/PIHALB,PIEINS,PI3HLB,PIZWEI
      K=J
      DO 1 I=1,N
      K=K+1
1     Y(K)=XX(I)
      DO 2 I=1,NS
      K=K+1
2     Y(K)=MAX(S(I),EPSILO(1))
      IF(.NOT.BKORRL) GOTO 4
      DO 3 I=1,NP
      K=K+1
      PI=P(I)
      IF(ABS(PI).GT.PIEINS) PI=PI-SIGN(PIZWEI,PI)
3     Y(K)=PI
4     K=K+1
      Y(K)=ZZ
      RETURN
      END
      SUBROUTINE UMSPEI
     1(K1,K2,KK,NY,Y)
      implicit none
      double precision Y(NY)
      integer K1,K2,KK,NY,K
      DO 1 K=1,KK
1     Y(K2+K)=Y(K1+K)
      RETURN
      END
      double precision FUNCTION Z(S,R)
      implicit none
      double precision ZP, A,B,D,S
      integer LZ
      double precision R
      external R
      COMMON/EVZ/LZ
      DATA ZP/6.28318531d0/
      if (LZ == 1) then
        A=SQRT(-2.* LOG(R(D)))
        B=ZP*R(D)
        Z=S*A*SIN(B)
        LZ=2
        RETURN
      else if (LZ == 2) then
        Z=S*A*COS(B)
        LZ=1
        RETURN
      else
        RETURN
      endif
      END
      double precision FUNCTION ZULASS
     1(N,M,XX,RESTRI)
      double precision XX(N)
      ZULASS=0.
      DO 1 J=1,M
      R=RESTRI(J,N,XX)
      IF(R.LT.0.) ZULASS=ZULASS-R
1     CONTINUE
      RETURN
      END
