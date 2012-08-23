SUBROUTINE rsigma(kanal,datei)

!!!$     Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

!!!$     Andreas Kemna                                            20-Dec-1993
!!!$     Letzte Aenderung   07-Nov-1997

!!!$.....................................................................

  USE make_noise
  USE alloci,ONLY:rnd_r,rnd_p
  USE femmod
  USE invmod
  USE sigmamod
  USE modelmod
  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT NONE

!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND=4) ::    kanal

!!!$     Datei
  CHARACTER (80)   ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     i,idum,ifp1,ifp2
  REAL(KIND(0D0))  ::     bet,pha,eps_r,eps_p
!!!$
  LOGICAL          ::     has_wmfak
!!!$ Pi
  REAL (KIND(0D0)) ::    pi
!!!$.....................................................................
  pi = dacos(-1d0)

  has_wmfak = .FALSE.
!!!$     'datei' oeffnen
  fetxt = datei

  errnr = 1
  OPEN(kanal,file=TRIM(fetxt),status='old',err=999)

  errnr = 3
  IF (lnsepri) THEN
     PRINT*,''
     CALL get_unit(ifp1)
     PRINT*,'iseedpri / std. err.',iseedpri,modl_stdn
     OPEN (ifp1,FILE='tmp.mynoise_mprior_rho',STATUS='replace')
     WRITE (*,'(A)',advance='no')' Initializing noise '
     ALLOCATE (rnd_r(elanz))
     CALL Random_Init(iseedpri)
     DO i=1,elanz
        rnd_r(i) = Random_Gauss()
     END DO
     IF (.NOT.ldc) THEN
        CALL get_unit(ifp2)
        OPEN (ifp2,FILE='tmp.mynoise_mprior_phase',STATUS='replace')
        ALLOCATE (rnd_p(elanz))
        CALL Random_Init(-iseedpri)
        DO i=1,elanz
           rnd_p(i) = Random_Gauss()
        END DO
     END IF
     PRINT*,''
  END IF

!!!$     Anzahl der Messwerte lesen
!!!$ also check if we may use individual errors or not
!!$  READ(kanal,*,END=1001,err=11) idum,has_wmfak
!!$  IF (has_wmfak) PRINT*,'---> Reference model regularization !'
!!$  GOTO 12
!!$11 PRINT*,' no reference model regularization '
!!$  BACKSPACE (kanal)

!!!$     Anzahl der Elemente (ohne Randelemente) einlesen
  READ(kanal,*,END=1001,err=1000) idum
!12 CONTINUE
!!!$     Ggf. Fehlermeldung
  IF (idum.NE.elanz) THEN
     fetxt = ' '
     errnr = 47
     GOTO 1000
  END IF
  
!!!$     Betrag und Phase (in mrad) des komplexen Widerstandes einlesen
  DO i=1,elanz
     
     IF (has_wmfak) THEN
        READ(kanal,*,END=1001,err=1002) bet,pha,wmfak(i)
     ELSE
        READ(kanal,*,END=1001,err=1000) bet,pha
     END IF

     pha      = 1d-3*pha

     IF (lprior) THEN 
        !     set prior model ...             
        IF (bet0 <= 0. .OR. &
             (.NOT.ldc.AND.dabs(pha0).GT.1d3*pi)) THEN
           fetxt = 'starting model incorrect '
           errnr = 91
           GOTO 999
        END IF

        IF (bet > 0d0) THEN  
           !     TODO: meaningful phase check.. 
           IF (lnsepri) THEN
              eps_r = 1d-2*modl_stdn * bet
              WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')rnd_r(i),eps_r,bet
              bet = bet + rnd_r(i) * eps_r
              WRITE(ifp1,'(G14.4)')bet
              IF (.NOT. ldc) THEN
                 eps_p = 1d-2*modl_stdn*dabs(pha)
                 WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')rnd_p(i),eps_p,pha
                 pha = pha + rnd_p(i) * eps_p
                 WRITE(ifp2,'(G14.4)')pha
              END IF
           END IF
           sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)
!!!$    hier koennte mittelung passieren
           m0(mnr(i)) = cdlog(sigma(i))
        ELSE                
           !     or let it stay at zero and take background cond
           sigma(i) = dcmplx( dcos(pha0/1d3)/bet0 , -dsin(pha0/1d3)/bet0 )
           m0(mnr(i)) = DCMPLX(0d0)
        END IF

     ELSE

!!!$     Ggf. Fehlermeldung
        IF (bet.LT.1d-12) THEN
           fetxt = ' '
           errnr = 11
           GOTO 999
        END IF

!!!$     Komplexe Leitfaehigkeit bestimmen
        sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)

     END IF
  END DO

!!!$     'datei' schliessen
  CLOSE(kanal)

  IF (lnsepri) THEN
     CLOSE (ifp1)
     IF (.NOT.ldc) CLOSE(ifp2)
  END IF

  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 RETURN

1000 CLOSE(kanal)
  RETURN

1001 CLOSE(kanal)
  errnr = 2
  RETURN

1002 CLOSE(kanal)
  fetxt = 'no reference factor'
  errnr = 2
  RETURN

END SUBROUTINE rsigma
