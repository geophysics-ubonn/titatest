subroutine rsigma(kanal,datei)

!!!$     Unterprogramm zum Einlesen der Widerstandsverteilung aus 'datei'.

!!!$     Andreas Kemna                                            20-Dec-1993
!!!$     Letzte Aenderung   07-Nov-1997

!!!$.....................................................................

  USE make_noise
  USE alloci,only:rnd_r,rnd_p
  USE femmod
  USE invmod
  USE sigmamod
  USE modelmod
  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT none

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

!!!$ Pi
  REAL (KIND(0D0)) ::    pi
!!!$.....................................................................
  pi = dacos(-1d0)

!!!$     'datei' oeffnen
  fetxt = datei

  errnr = 1
  open(kanal,file=fetxt,status='old',err=999)

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

!!!$     Anzahl der Elemente (ohne Randelemente) einlesen
  read(kanal,*,end=1001,err=1000) idum

!!!$     Ggf. Fehlermeldung
  if (idum.ne.elanz) then
     fetxt = ' '
     errnr = 47
     goto 1000
  end if

!!!$     Betrag und Phase (in mrad) des komplexen Widerstandes einlesen
  do i=1,elanz
     read(kanal,*,end=1001,err=1000) bet,pha

     pha      = 1d-3*pha

     IF (lprior) THEN 
        !     set prior model ...             
        IF (bet0 <= 0. .OR. &
             (.not.ldc.and.dabs(pha0).gt.1d3*pi)) THEN
           fetxt = 'starting model incorrect '
           errnr = 91
           goto 999
        END IF

        IF (bet > 0.) THEN  
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
        if (bet.lt.1d-12) then
           fetxt = ' '
           errnr = 11
           goto 999
        end if

!!!$     Komplexe Leitfaehigkeit bestimmen
        sigma(i) = dcmplx(dcos(pha)/bet,-dsin(pha)/bet)

     END IF
  end do

!!!$     'datei' schliessen
  close(kanal)

  IF (lnsepri) THEN
     close (ifp1)
     IF (.NOT.ldc) close(ifp2)
  END IF

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

1001 close(kanal)
  errnr = 2
  return

end subroutine rsigma
