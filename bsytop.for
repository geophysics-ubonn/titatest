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

      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin'        ! fuer sytop und den ganzen rest 

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

      iel = 0
      DO i=1,typanz

         iel = iel + nelanz(i)

         IF (typ(i) /= 11) CYCLE  ! suche nach "no flow"

         nkel = selanz(i)
         sytop = 0.

         DO j=1,nelanz(i)

            iel = iel + 1

            sp=0
            DO ik=1,nkel
               sp = sp + sy(snr(nrel(iel,ik)))
            END DO
            sp = sp/DBLE(nkel)


            IF (j == 1) THEN !errechne direkt den mittelwert
               sytop = sp
            ELSE
               sytop = (sytop*dble(j-1)+sp)/dble(j)
            END IF

         END DO
      END DO

      END SUBROUTINE bsytop
      
