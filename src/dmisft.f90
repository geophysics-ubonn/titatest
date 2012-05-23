SUBROUTINE dmisft(lsetup)

!!!$     Unterprogramm zum Bestimmen des Misfits der Daten.

!!!$     Andreas Kemna                                            01-Mar-1995
!!!$     Letzte Aenderung   15-Jan-2001

!!!$.....................................................................

  USE invmod
  USE datmod
  USE errmod
  USE konvmod
  USE pathmod

  IMPLICIT NONE


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Hilfsschalter
  LOGICAL ::     lsetup

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfelder
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: psi,eps2
  INTEGER(KIND = 4),DIMENSION(:),ALLOCATABLE :: wdlok

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     i,idum
  COMPLEX (KIND(0D0)) ::    cdum,cdat,csig
  REAL(KIND(0D0))     ::     dum,norm,norm2

!!!$.....................................................................

!!!$     einfach mal oeffnen falls Ausgabe
  errnr = 1
  fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
  OPEN(fpeps,file=TRIM(fetxt),STATUS='old',POSITION='append',ERR=1000)
  fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
  OPEN(fprun,file=TRIM(fetxt),STATUS='old',POSITION='append',ERR=1000)
  errnr = 4

11 FORMAT (I6,G12.3,1x,F10.6,2x,I3,2X,4(1X,G15.7))

  IF ((llam.AND..NOT.lstep).OR.lsetup) THEN
     IF (lip) THEN
        WRITE(fpeps,'(/a,t7,I4)',err=1000)'PIT#',it
     ELSE
        WRITE(fpeps,'(/a,t7,I4)',err=1000)'IT#',it
     END IF
     WRITE(fpeps,'(a)',err=1000)'Datum'//ACHAR(9)//'eps'//ACHAR(9)//&
          'psi'//ACHAR(9)//'pol'//ACHAR(9)//'Re(d)'//ACHAR(9)//'Re(f(m))'&
          //ACHAR(9)//'Im(d)'//ACHAR(9)//'Im(f(m))'
  END IF

!!!$     RMS-WERTE BERECHNEN
  nrmsd  = 0d0
  betrms = 0d0
  pharms = 0d0
  idum   = 0
  !     get memory for wdlok and psi
  ALLOCATE (wdlok(nanz),psi(nanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation psi'
     errnr = 94
     RETURN
  END IF

  DO i=1,nanz
     wdlok(i) = 1

!!!$     Phasen lokal korrigieren
     CALL chkpo2(dat(i),sigmaa(i),cdat,csig,wdlok(i),lpol)

!!!$     ak Ggf. Daten mit Phase betraglich groesser 200 mrad nicht beruecksichtigen
!!!$     ak (Standardabweichung auf 1d4 hochsetzen)
!!!$     ak           if (dabs(1d3*dimag(csig)).gt.200d0) wmatd(i)=1d-8

!!!$     diff-            cdum = cdat - csig
!!!$     diff+<
     IF (.NOT.ldiff) THEN
        cdum = cdat - csig
     ELSE
        cdum = cdat - csig - (d0(i) - fm0(i))
     END IF
!!!$     diff+>

     IF (lip) THEN
        psi(i) = dsqrt(wmatd(i))*dabs(dimag(cdum))
     ELSE
        psi(i) = dsqrt(wmatd(i))*cdabs(cdum)
     END IF

!!!$     Ggf. 'eps_i', 'psi_i' und Hilfsfeld ausgeben
     IF ((llam.AND..NOT.lstep).OR.lsetup) WRITE(fpeps,11,err=1000) &
          i, REAL(1d0/dsqrt(wmatd(i))),REAL(psi(i)),wdlok(i),&
          REAL(csig),REAL(cdat),AIMAG(cdat),AIMAG(csig)

     idum   = idum   + wdlok(i)
     nrmsd  = nrmsd  + psi(i)*psi(i)*DBLE(wdlok(i))

     betrms = betrms + wmatdr(i)*DBLE(wdlok(i)) * &
          DBLE(cdum)*DBLE(cdum)

     pharms = pharms + wmatdp(i)*DBLE(wdlok(i)) * &
          dimag(cdum)*dimag(cdum)

  END DO

!!!$     Ggf. Fehlermeldung
  IF (idum.EQ.0) THEN
     fetxt = ' '
     errnr = 99
     GOTO 1000
  END IF

  npol   = nanz-idum
  rmssum = nrmsd
  nrmsd  = dsqrt(nrmsd /DBLE(idum))
  betrms = dsqrt(betrms/DBLE(idum))
  pharms = dsqrt(pharms/DBLE(idum))

!!!$     Ggf. ROBUST INVERSION (nach Doug' LaBrecque)
  IF (lrobust) THEN

     !     get memory for wdlok
     ALLOCATE (eps2(nanz),stat=errnr)
     IF (errnr /= 0) THEN
        fetxt = 'Error memory allocation eps2'
        errnr = 94
        RETURN
     END IF
!!!$     'estimated weights' und 1-Normen berechnen
     norm  = 0d0
     norm2 = 0d0

     DO i=1,nanz
        dum     = 1d0/SQRT(wmatd(i))
        eps2(i) = dum*SQRT(psi(i))
        norm    = norm  + psi(i)*(wdlok(i))
        norm2   = norm2 + psi(i)*(wdlok(i))*dum/eps2(i)
     END DO

!!!$     'estimated weights' normieren
     DO i=1,nanz
        eps2(i) = eps2(i) * norm2/norm
     END DO

!!!$     Kleinere Standardabweichung ausschliessen und neue 1-Norm berechnen
     norm2 = 0d0

     DO i=1,nanz
        dum     = 1d0/SQRT(wmatd(i))
        eps2(i) = MAX(dum,eps2(i))
        norm2   = norm2 + (psi(i))*(wdlok(i))*dum/eps2(i)
     END DO

     l1rat = norm/norm2

!!!$     Ggf. neue Wichtungsfaktoren belegen
     IF (l1rat.GT.l1min) THEN
        DO i=1,nanz
           dum = 1./(eps2(i) * eps2(i))

!!!$     Ausgabe, falls 'eps_neu' > 1.1 * 'eps_alt'
           IF (dum.LT.0.83d0*wmatd(i).AND. &
                ((llam.AND..NOT.lstep).OR.lsetup)) THEN

              WRITE(fprun,*)i,&
                   ' : increase standard deviation *', &
                   REAL(eps2(i)*SQRT(wmatd(i))),dum,psi(i)
           END IF

           wmatd(i) = dum
        END DO
     END IF
     IF (ALLOCATED(eps2)) DEALLOCATE (eps2)
  END IF

  IF (ALLOCATED(wdlok)) DEALLOCATE (wdlok,psi)

  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen
1000 CLOSE(fpeps)
  CLOSE(fprun)
  RETURN

END SUBROUTINE dmisft


!!!$*********************************************************************

SUBROUTINE chkpo2(dati,sigi,cdat,csig,wdlok,ldum)

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Eingabe
  COMPLEX (KIND(0D0)) ::    dati,sigi

!!!$     Ausgabe
  COMPLEX (KIND(0D0)) ::    cdat,csig

!!!$     Schalter
  INTEGER (KIND = 4)  ::     wdlok
  LOGICAL ::     ldum

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND = 4)  ::     idat,isig

!!!$     Real-, Imaginaerteile
  REAL (KIND(0D0))    ::     redat,imdat,resig,imsig

!!!$     Pi
  REAL (KIND(0D0))    ::     pi

!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     Logarithmierte Betraege in den Realteilen,
!!!$     Phasen (in rad) in den Imaginaerteilen
  redat = DBLE(dati)
  imdat = dimag(dati)
  resig = DBLE(sigi)
  imsig = dimag(sigi)

!!!$     Phasenbereich checken
  IF (imdat.GT.pi/2d0) THEN
     idat = -1
  ELSE IF (imdat.LE.-pi/2d0) THEN
     idat = 1
  ELSE
     idat = 0
  END IF

  IF (imsig.GT.pi/2d0) THEN
     isig = -1
  ELSE IF (imsig.LE.-pi/2d0) THEN
     isig = 1
  ELSE
     isig = 0
  END IF

  IF (idat.EQ.0.AND.isig.NE.0) THEN

!!!$     Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
!!!$     ggf. Korrektur; auf jeden Fall Polaritaetswechsel
     imsig = imsig + DBLE(isig)*pi
     IF (.NOT.ldum) imdat=imdat-dsign(pi,imdat)

     wdlok = 0

  ELSE IF (idat.NE.0.AND.isig.EQ.0) THEN

!!!$     Falls ldum=.true., angenommene Polaritaet des Messdatums falsch,
!!!$     ggf. Korrektur
     IF (ldum) imdat=imdat+DBLE(idat)*pi

     wdlok = 0

  ELSE IF (idat.NE.0.AND.isig.NE.0) THEN

!!!$     Polaritaetswechsel
     imsig = imsig + DBLE(isig)*pi
     imdat = imdat + DBLE(idat)*pi

  END IF

!!!$     'cdat' und 'csig' speichern
  cdat = dcmplx(redat,imdat)
  csig = dcmplx(resig,imsig)

  RETURN
END SUBROUTINE chkpo2
