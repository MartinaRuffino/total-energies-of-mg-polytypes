      SUBROUTINE SPCFND(SID,SPECID,NSPEC,J)
      CHARACTER*(*) SID,SPECID(1)
      J=0
      DO 10 I=1,NSPEC
      IF(SID.EQ.SPECID(I)) THEN
        J=I
        RETURN
        ENDIF
  10  CONTINUE
      RETURN
      END
