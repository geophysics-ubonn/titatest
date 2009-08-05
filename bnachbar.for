      SUBROUTINE bnachbar
c     
c     Unterprogramm zum Bestimmen der Elementnachbarn
c     zur Realisierung der Triangulationsregularisierung
c     in CRTomo von Andreas Kemna
c     
c     Roland Martin                                            29-Jul-2009
c     
c     Letzte Aenderung                                         29-Jul-2009
c.....................................................................
      

      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin'        ! fuer nachbar, nrel etc. 
      INCLUDE 'model.fin'       ! fuer manz

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariablen
      INTEGER :: i,j,ik,jk
c     Knotennummer der Kanten zaehler von element i und j
      INTEGER :: ik1,ik2,jk1,jk2
c     maximale knotenanzahl nicht entarteter elemente
      INTEGER :: smaxs
c---------------------------------------------
      
      smaxs = selanz(1)

      IF (smaxs /= MAXVAL(selanz)) THEN
         PRINT*,'smaxs/=MAXVAL(selanz)!! check grid file please'
         STOP
      END IF

      nachbar = 0

      DO i=1,manz

         DO ik=1,smaxs

            ik1 = nrel(i,ik)
            ik2 = nrel(i,MOD(ik,smaxs)+1)

            DO j=1,manz

               IF (j==i) CYCLE
               
               DO jk=1,smaxs
                  
                  jk1 = nrel(j,jk)
                  jk2 = nrel(j,MOD(jk,smaxs)+1)
                  
                  IF ( (ik1==jk1.AND.ik2==jk2) .OR.
     1                 (ik1==jk2.AND.ik2==jk1) ) THEN
                     
                     nachbar(i,ik) = j ! Element teilt kante
                     nachbar(i,0) = nachbar(i,0)+1 ! Anzahl der Nachbarn
                     
                  END IF
                  
               END DO
            END DO              ! inner loop j=1,manz

         END DO
      END DO                    ! outer loop i=1,manz


      END
