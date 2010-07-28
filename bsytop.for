      SUBROUTINE bsytop
c     
c     Unterprogramm zum Bestimmen der Mittleren Elementhoehe
c     aller "no-flow" Elemente (Halbraumgrenze) sytop
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          20-Nov-2009
c     
c     Letzte Aenderung                                         20-Nov-2009
c     
c.....................................................................

      USE elemmod        ! fuer sytop und den ganzen rest 

      IMPLICIT none

c     PROGRAMMINTERNE PARAMETER:-------------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik
c     Elementnummer
      INTEGER :: iel
c     Knotenanzahl der no flow elemente 
      INTEGER :: nkel
c     Schwerpunktskoordinaten des randelements
      REAL(KIND(0D0)) :: sp
c-----------------------------------------------------------------------

      sytop = 0.
      iel = 0
      DO i=1,typanz

         iel = iel + nelanz(i)

         IF (typ(i) /= 12) CYCLE  ! suche nach "no flow"

         nkel = selanz(i)
         sytop = 0.

         DO j=1,nelanz(i)

            iel = iel + 1

            sp=0
            DO ik=1,nkel
               PRINT*,ik,sp,sy(snr(nrel(iel,ik))),sx(snr(nrel(iel,ik)))
               sp = sp + sy(snr(nrel(iel,ik)))
            END DO
            sp = sp/DBLE(nkel)


            IF (j == 1) THEN 
               sytop = sp
            ELSE                !errechne den mittelwert (diret)
               sytop = (sytop*dble(j-1)+sp)/dble(j)
            END IF

!            PRINT*,sp,sytop,j

         END DO
      END DO

!      PRINT*,sytop


      END SUBROUTINE bsytop
      
