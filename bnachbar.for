      SUBROUTINE bnachbar
c     
c     Unterprogramm zum Bestimmen der Elementnachbarn
c     zur Realisierung der Triangulationsregularisierung
c     sowie der kleinsten moeglichen Skalenlaenge der
c     stochastischen Regularisierung
c     
c     Copyright by Andreas Kemna 2009
c     
c     Erste Version von Roland Martin                          29-Jul-2009
c     
c     Letzte Aenderung                                         07-Aug-2009
c     
c.....................................................................

      USE alloci
      USE modelmod       ! fuer manz
      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin'        ! fuer nachbar, nrel etc. 
      INCLUDE 'err.fin'       ! fuer manz

c     PROGRAMMINTERNE PARAMETER:----------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik,jk
c     Knotennummer der Kanten zaehler von element i und j
      INTEGER :: ik1,ik2,jk1,jk2
c     Maximale knotenanzahl nicht entarteter Elemente
      INTEGER :: smaxs
c-----------------------------------------------------------------------

      smaxs = MAXVAL(selanz)

      IF (.NOT.ALLOCATED (nachbar)) 
     1     ALLOCATE (nachbar(manz,smaxs+1),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem nachbar in bnachbar'
         errnr = 97
         RETURN
      END IF

      nachbar = 0

      DO i=1,elanz
         
         WRITE (*,'(a,1X,F6.2,a)',ADVANCE='no')ACHAR(13)//
     1        'bnachbar/ ',REAL (i * (100./elanz)),'%'

         DO ik=1,smaxs

            ik1 = nrel(i,ik)
            ik2 = nrel(i,MOD(ik,smaxs)+1)

            DO j=1,elanz

               IF (j==i) CYCLE

               DO jk=1,smaxs

                  jk1 = nrel(j,jk)
                  jk2 = nrel(j,MOD(jk,smaxs)+1)

                  IF ( (ik1==jk1.AND.ik2==jk2) .OR.
     1                 (ik1==jk2.AND.ik2==jk1) ) THEN

                     nachbar(i,ik) = j ! Element teilt kante
                     nachbar(i,smaxs+1) = nachbar(i,smaxs+1)+1 ! Anzahl der Nachbarn

                  END IF

               END DO
            END DO              ! inner loop j=1,elanz

         END DO
      END DO                    ! outer loop i=1,elanz

      END SUBROUTINE bnachbar
