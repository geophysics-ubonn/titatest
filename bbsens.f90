subroutine bbsens(kanal,datei)

!!!$     Unterprogramm zur Berechnung der Summe der Sensitivitaeten 
!!!$     aller Messungen (normiert)
!!!$     Berechnet nun coverage als Summe der Absolutbetraege..
!!!$     ('kpot' als Hilfsfeld benutzt).

!!!$     Andreas Kemna                                   02-Mar-1995
!!!$     Letzte Aenderung              31-Mar-2010

!!!$....................................................................

  USE alloci , ONLY : sens,sensdc
  USE datmod , ONLY : nanz
  USE invmod , ONLY : lip,wmatd,wdfak
  USE modelmod , ONLY : manz
  USE elemmod , ONLY: espx,espy
  USE errmod , ONLY : errnr,fetxt
  USE femmod , ONLY : ldc
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
  REAL (KIND(0D0)),ALLOCATABLE,DIMENSION(:) :: csens
  REAL (KIND(0D0))                          :: csensmax
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

  IF (lip) THEN
     DO j=1,manz
        DO i=1,nanz
           csens(j) = csens(j) + DBLE(sens(i,j)) * &
                DBLE(sens(i,j)) * wmatd(i)*DBLE(wdfak(i))
        END DO
     END DO
  ELSE IF (ldc) THEN
     DO j=1,manz
        DO i=1,nanz
           csens(j) = csens(j) + sensdc(i,j) * &
                sensdc(i,j) * wmatd(i)*DBLE(wdfak(i))
        END DO
     END DO
  ELSE
     DO j=1,manz
        DO i=1,nanz
           csens(j) = csens(j) + DCONJG(sens(i,j)) * &
                sens(i,j) * wmatd(i)*dble(wdfak(i)) 
!!!$ wechselt automatisch wmatdp bei lip
        END DO
     END DO
  ENDIF

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
  
  DEALLOCATE (csens)
  
  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

end subroutine bbsens
