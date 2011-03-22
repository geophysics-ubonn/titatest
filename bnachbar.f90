SUBROUTINE bnachbar
!!!$     
!!!$     Unterprogramm zum Bestimmen der Elementnachbarn
!!!$     zur Realisierung der Triangulationsregularisierung
!!!$     sowie der kleinsten moeglichen Skalenlaenge der
!!!$     stochastischen Regularisierung
!!!$     
!!!$     Copyright by Andreas Kemna 2009
!!!$     
!!!$     Erste Version von Roland Martin                          29-Jul-2009
!!!$     
!!!$     Letzte Aenderung                                         07-Aug-2009
!!!$     
!!!$.....................................................................

  USE alloci
  USE modelmod       ! fuer manz
  USE elemmod        ! fuer nachbar, nrel etc. 
  USE errmod       ! errnr und fetxt
  USE konvmod , ONLY : lverb

  IMPLICIT none


!!!$     PROGRAMMINTERNE PARAMETER:----------------------------------------
!!!$     Indexvariablen
  INTEGER :: i,j,ik,jk,count
!!!$     Knotennummer der Kanten zaehler von element i und j
  INTEGER :: ik1,ik2,jk1,jk2
!!!$-----------------------------------------------------------------------

  IF (.NOT.ALLOCATED (nachbar)) &
       ALLOCATE (nachbar(manz,smaxs+1),STAT=errnr)
  IF (errnr/=0) THEN
     WRITE (*,'(/a/)')'Allocation problem nachbar in bnachbar'
     errnr = 97
     RETURN
  END IF

  nachbar = 0
  count = 0
  !$OMP PARALLEL DEFAULT (none) &
  !$OMP PRIVATE (i,ik,ik1,ik2,j,jk,jk1,jk2) &
  !$OMP SHARED (nachbar,smaxs,count,elanz,nrel,lverb)
  !$OMP DO
  DO i=1,elanz

     !$OMP ATOMIC
     count = count + 1

     IF (lverb) WRITE (*,'(a,t70,F6.2,a)',ADVANCE='no')ACHAR(13)// &
          'bnachbar/ ',REAL (count * (100./elanz)),'%'

     DO ik=1,smaxs

        ik1 = nrel(i,ik)
        ik2 = nrel(i,MOD(ik,smaxs)+1)

        DO j=1,elanz

           IF (j==i) CYCLE

           DO jk=1,smaxs

              jk1 = nrel(j,jk)
              jk2 = nrel(j,MOD(jk,smaxs)+1)

              IF ( (ik1==jk1.AND.ik2==jk2) .OR. &
                   (ik1==jk2.AND.ik2==jk1) ) THEN

                 nachbar(i,ik) = j ! Element teilt kante
                 nachbar(i,smaxs+1) = nachbar(i,smaxs+1)+1 ! Anzahl der Nachbarn

              END IF

           END DO
        END DO              ! inner loop j=1,elanz

     END DO
  END DO                    ! outer loop i=1,elanz
  !$OMP END PARALLEL

END SUBROUTINE bnachbar
