      SUBROUTINE MDIAN1(X,N,XMED)
      IMPLICIT none
      REAL(KIND(0D0)),DIMENSION(N) :: X
      REAL(KIND(0D0)) :: XMED
      INTEGER :: N,N2

      CALL SORT(N,X)
      N2=N/2
      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
      ELSE
        XMED=X(N2+1)
      ENDIF
      RETURN
      END
