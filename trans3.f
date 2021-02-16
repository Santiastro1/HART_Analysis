C ______ Transform real X.XXXX to character*6______


      SUBROUTINE TRANS3(CI,XN)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*5 CI
      real*4 xn 
c a los primeros 5 caracteres le escribe el real, el real tiene un
c decimal.
      WRITE(CI,'(F5.3)')XN



c al resto de les caracteres los ignora.
      DO I=1,4
         IF(CI(I:I).EQ.' ')THEN
           CI(I:I)='0'
         ENDIF
      ENDDO
      RETURN
      END

