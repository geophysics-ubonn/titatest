subroutine bbsens(kanal,datei)

!!!$     Unterprogramm zur Berechnung der Summe der Sensitivitaeten 
!!!$     aller Messungen (normiert)
!!!$     Berechnet nun coverage als Summe der Absolutbetraege..
!!!$     ('kpot' als Hilfsfeld benutzt).

!!!$     Andreas Kemna                                   02-Mar-1995
!!!$     Letzte Aenderung              31-Mar-2010

!!!$....................................................................

  USE alloci , ONLY : sens,sensdc, csens
  USE datmod , ONLY : nanz
  USE invmod , ONLY : lip,wmatd,wdfak
  USE modelmod , ONLY : manz
  USE elemmod , ONLY: espx,espy
  USE errmod , ONLY : errnr,fetxt
  USE femmod , ONLY : ldc
  USE ompmod
  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4)  ::     kanal

!!!$     Datei
  CHARACTER (80)      ::  datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
!!$  REAL (KIND(0D0)),ALLOCATABLE,DIMENSION(:) :: csens
  REAL (KIND(0D0))                          :: csensmax
  REAL (KIND(0D0))                          :: dum
!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j
!!!$.....................................................................


!!!$     'datei' oeffnen
  fetxt = datei
  errnr = 1
  open(kanal,file=TRIM(fetxt),status='replace',err=999)
  errnr = 4
  ALLOCATE (csens(manz),STAT=errnr)
  IF (errnr /= 0) RETURN
!!!$     Werte berechnen
  csens = 0D0
  !$OMP PARALLEL &
  !$OMP SHARED (manz,csens,wdfak,wmatd,sens,sensdc,lip,ldc,nanz) &
  !$OMP PRIVATE (dum,i)
  !$OMP DO SCHEDULE (GUIDED,CHUNK_0)
  DO j=1,manz
     DO i=1,nanz
        dum = SQRT(wmatd(i)) * DBLE(wdfak(i))
        IF (lip) THEN
           csens(j) = csens(j) + ABS(DBLE(sens(i,j))) * dum
        ELSE IF (ldc) THEN
           csens(j) = csens(j) + ABS(sensdc(i,j)) * dum
        ELSE
!!!$ wechselt automatisch wmatdp bei lip
           csens(j) = csens(j) + ABS(sens(i,j)) * dum
        ENDIF
     END DO
  END DO
  !$OMP END PARALLEL
!!!$ for normalization
  csensmax = MAXVAL(csens)
  
  write(kanal,*,err=1000) manz

!!!$     Koordinaten und Sensitivitaetsbetraege schreiben
!!!$     (logarithmierter (Basis 10) normierter Betrag)
  DO i=1,manz
     WRITE (kanal,*,err=1000)espx(i),espy(i),LOG10(csens(i)/csensmax)
  END DO

!!!$     Maximale Sensitivitaet schreiben
  write(kanal,*,err=1000)'Max:',csensmax
!!!$     'datei' schliessen
  close(kanal)
  
!!$  DEALLOCATE (csens)
  
  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

end subroutine bbsens
