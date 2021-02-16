c ------------------------------------------------


      SUBROUTINE TRANS2(CI,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*3 CI
   
      WRITE(CI,'(I3)')N

      DO I=1,3
         IF(CI(I:I).EQ.' ')THEN
           CI(I:I)='0'
         ENDIF
      ENDDO
      RETURN
      END

