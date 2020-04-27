      SUBROUTINE IUNPAC(T,L,NL)
C  UNPACKS STRING OF INTEGER VALUES
      DIMENSION L(10)
      CHARACTER*(*) T
      CHARACTER*1 DIG(-1:9)
      DATA DIG/'-','0','1','2','3','4','5','6','7','8','9'/
      NT=LEN(T)
      NL=0
      DO 10 J=1,NT
      IF(T(J:J).EQ.' ') RETURN
      DO 11 I=-1,9
      IF(T(J:J).EQ.DIG(I)) THEN
         NL=NL+1
         L(NL)=I
         GOTO 10
         ENDIF
  11  CONTINUE
      WRITE(6,*) '----  IUNPAC:  INPUT=<',T,'>'
      CALL RX('STRING HAS WRONG FORMAT''')
  10  CONTINUE
      RETURN
      END
