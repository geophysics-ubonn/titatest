SUBROUTINE bsytop
!!!$     
!!!$     Unterprogramm zum Bestimmen der Mittleren Elementhoehe
!!!$     aller "no-flow" Elemente (Halbraumgrenze) sytop
!!!$     
!!!$     Copyright by Andreas Kemna
!!!$     
!!!$     Erste Version von Roland Martin                          20-Nov-2009
!!!$     
!!!$     Letzte Aenderung                                         20-Nov-2009
!!!$     
!!!$.....................................................................
use alloci, only:prec
  USE elemmod               ! fuer sytop und den ganzen rest 

  IMPLICIT none

!!!$     PROGRAMMINTERNE PARAMETER:-------------------------------------------
!!!$     Indexvariablen
  INTEGER :: i,j,ik
!!!$     Elementnummer
  INTEGER :: iel
!!!$     Knotenanzahl der no flow elemente 
  INTEGER :: nkel
!!!$     Schwerpunktskoordinaten des randelements
  REAL(prec) :: sp
!!!$-----------------------------------------------------------------------

  sytop = 0.
  iel = 0
  DO i=1,typanz

     IF (typ(i) == 12) THEN ! search for "no flow" boundary

        nkel = selanz(i)
        sytop = 0.

        DO j=1,nelanz(i) ! check all electrode y-positions

           iel = iel + 1

           sp=0
           DO ik=1,nkel
              sp = sp + sy(snr(nrel(iel,ik)))
           END DO
           sp = sp/REAL(nkel)


           IF (j == 1) THEN 
              sytop = sp
           ELSE             !caculate arithmetic mean (direct)
              sytop = (sytop*REAL(j-1)+sp)/REAL(j)
           END IF

        END DO
     ELSE

        iel = iel + nelanz(i)

     END IF
  END DO

  IF (sytop /= 0.) PRINT*,'SYTOP::',sytop

END SUBROUTINE bsytop

