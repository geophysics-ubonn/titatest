      SUBROUTINE bnachbar
      
c     Unterprogramm zum Bestimmen der Elementnachbarn
c     zur Realisierung der Triangulationsregularisierung
c     in CRTomo.
c
c     Andreas Kemna
ci.A. Roland Martin                                            29-Jul-2009
c     Letzte Aenderung   29-Jul-2009
c 
c.....................................................................
 

      IMPLICIT none

      INCLUDE 'parmax.fin' ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin' ! fuer nachbar, nrel etc. 

c PROGRAMMINTERNE PARAMETER:

c Indexvariablen
      integer         * 4     i,j,k,ik,j2
c Knotennummer von Kanten zähler
      integer         * 4     ik1,ik2,in1,in2

      nachbar=0

      DO i=1,typanz

         IF (typ(i)>10) EXIT    ! nur die flächenelemente die als erste kommen

         DO j=1,nelanz(i)       ! alle FL-elemente durchgehen

            nachbar(j,0)=selanz(i) ! knotenpunkte speichern :(

            DO k=1,selanz(i)    ! alle kanten jeden Elements

               ik1=nrel(j,k)
               IF (k==selanz(i)) THEN
                  ik2=nrel(j,1)
               ELSE
                  ik2=nrel(j,k+1)
               END IF

               DO j2=1,nelanz(i) ! suche nach gemeinsamen Knoten

                  IF (j2==j) CYCLE

                  DO ik=1,selanz(i)
                     
                     in1=nrel(j2,ik)

                     IF (ik==selanz(i)) THEN
                        in2=nrel(j2,1)
                     ELSE
                        in2=nrel(j2,ik+1)
                     END IF

                     IF ( (ik1==in1.AND.ik2==in2).OR.
     1                    (ik1==in2.AND.ik2==in1)) nachbar(j,k)=j2

                  END DO
               END DO
            END DO
         END DO
      END DO

      END 
