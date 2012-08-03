SUBROUTINE rdati(kanal,datei)

!!!$     Unterprogramm zum Einlesen der Elektrodenkennungen und der Daten
!!!$     inkl. Standardabweichungen aus 'datei'.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   20-Aug-2007

!!!$.....................................................................
  USE make_noise
  USE alloci, ONLY:rnd_r,rnd_p
  USE femmod
  USE datmod
  USE invmod
  USE electrmod
  USE errmod
  USE konvmod

  IMPLICIT NONE


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4) ::     kanal

!!!$     Datei
  CHARACTER (80)     ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND = 4) ::     i,ifp1,ifp2,ifp3

!!! >> RM
!!!$  USED ONLY IF NOISE IS ADDED!!
!!!$ Magnitude and Phase of new data
  REAL(KIND(0D0))     ::     new_bet,new_pha
!!!$ Counting signum mismatches of magnitude and phases 
!!!$    if noise is added
  INTEGER (KIND = 4) ::     icount_pha,icount_mag
!!!$ << RM

!!!$     Elektrodennummern
  INTEGER (KIND = 4) ::     elec1,elec2,elec3,elec4

!!!$     Betrag und Phase (in mrad) der Daten
  REAL(KIND(0D0))     ::     bet,pha

!!!$     Standardabweichung eines logarithmierten (!) Datums
  REAL(KIND(0D0))      ::     stabw
!!!$     Error of the resistance
  REAL(KIND(0D0))     ::     eps_r
!!!$     Standardabweichung der Phase
  REAL(KIND(0D0))     ::     stabwp,stabwb
!!!$     Error of the phase
  REAL(KIND(0D0))     ::     eps_p
!!!$     Pi
  REAL(KIND(0D0))     ::     pi
!!!$ check whether the file format is crtomo konform or not..
  LOGICAL             ::    crtf
!!!$.....................................................................
  pi = dacos(-1d0)

!!!$     'datei' oeffnen
  fetxt = datei
  errnr = 1
  OPEN(kanal,file=TRIM(fetxt),status='old',err=999)
  errnr = 3

!!!$     Anzahl der Messwerte lesen
!!!$ also check if we may use individual errors or not
  READ(kanal,*,END=1001,err=11) nanz,lindiv
  IF (lindiv) PRINT*,'+ Individual data error!'
  GOTO 12
11 PRINT*,'+ Taking error model'
  BACKSPACE(kanal)
!!!$c check if data file format is CRTOmo konform..
12 READ(kanal,*,END=1001,err=1000) elec1
  BACKSPACE(kanal)

  elec3=elec1-10000 ! are we still positive?
  crtf=(elec3 > 0) ! crtomo konform?

  ALLOCATE (strnr(nanz),strom(nanz),volt(nanz),sigmaa(nanz),&
       kfak(nanz),wmatdr(nanz),wmatdp(nanz),vnr(nanz),dat(nanz),&
       wmatd(nanz),wmatd2(nanz),sgmaa2(nanz),wdfak(nanz),&
       stat=errnr)
  wmatd = 0.;wmatdp = 0.; wmatdr = 0.
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation data space'
     errnr = 94
     GOTO 1000
  END IF

!!!$     Stromelektrodennummern, Spannungselektrodennummern, Daten inkl.
!!!$     auf 1 normierte Standardabweichungen lesen und Daten logarithmieren
  IF (lnse ) THEN
     icount_mag = 0;icount_pha = 0
     WRITE (*,'(A)',ADVANCE='no')ACHAR(13)//'Initializing noise'
     CALL get_unit(ifp3)
     OPEN (ifp3,FILE='inv.mynoise_voltages',STATUS='replace')
     WRITE (ifp3,*) nanz
     CALL get_unit(ifp1)
     OPEN (ifp1,FILE='inv.mynoise_rho',STATUS='replace')
     WRITE(ifp1,'(a)')'#  rnd_r'//ACHAR(9)//'eps_r'//&
          ACHAR(9)//ACHAR(9)//'bet(old)'//ACHAR(9)//'bet(new)'
     IF (.NOT. ldc) THEN
        CALL get_unit(ifp2)
        OPEN (ifp2,FILE='inv.mynoise_pha',STATUS='replace')
        WRITE(ifp2,'(a)')'#  rnd_p'//ACHAR(9)//'eps_p'//&
             ACHAR(9)//ACHAR(9)//'pha(old)'//ACHAR(9)//'pha(new)'
     END IF

     ALLOCATE (rnd_r(nanz))
     CALL Random_Init(iseed)
     DO i=1,nanz
        rnd_r(i) = Random_Gauss()
     END DO
     IF (.NOT.ldc) THEN
        ALLOCATE (rnd_p(nanz))
        CALL Random_Init(-iseed)
        DO i=1,nanz
           rnd_p(i) = Random_Gauss()
        END DO
     END IF
  END IF


  DO i=1,nanz
     stabwp = 0.; stabwb = 0.
     IF (lverb) WRITE (*,'(A,t70,F6.2,A)',ADVANCE='no')ACHAR(13)//&
          'data set ',REAL( i * (100./nanz) ),'%'
     IF (lindiv) THEN
        IF (ldc) THEN
           IF (crtf) THEN
              READ(kanal,*,END=1001,err=1000)strnr(i),vnr(i),bet,stabw
           ELSE
              READ(kanal,*,END=1001,err=1000)elec1,elec2,elec3,elec4,bet,&
                   stabw
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF
        ELSE
           IF (lfphai) THEN
              IF (crtf) THEN
                 READ(kanal,*,END=1001,err=1000)strnr(i),vnr(i),bet,pha,&
                      stabw,stabwp
              ELSE
                 READ(kanal,*,END=1001,err=1000)elec1,elec2,elec3,elec4,&
                      bet,pha,stabw,stabwp
                 strnr(i) = elec1*10000 + elec2
                 vnr(i)   = elec3*10000 + elec4
              END IF
!!!$     Ggf. Fehlermeldung
              IF (stabwp.LE.0d0) THEN
                 fetxt = ' '
                 errnr = 88
                 GOTO 1000
              END IF

!!!$     ak                        stabwp = stabp0 * stabwp
              stabwp = 1d-3*stabwp
           ELSE
              IF (crtf) THEN
                 READ(kanal,*,END=1001,err=1000)strnr(i),vnr(i),bet,pha,&
                      stabw
              ELSE
                 READ(kanal,*,END=1001,err=1000)elec1,elec2,elec3,elec4,&
                      bet,pha,stabw
                 strnr(i) = elec1*10000 + elec2
                 vnr(i)   = elec3*10000 + elec4

              END IF
           END IF
        END IF

!!!$     Ggf. Fehlermeldung
        IF (stabw.LE.0d0) THEN
           fetxt = ' '
           errnr = 88
           GOTO 1000
        END IF

!!!$     ak                stabw = stabw0 * stabw

!!!$     ak!!!$         Ggf. Fehlermeldung
!!!$     ak                if (bet.le.0d0) then
!!!$     akcak
!!!$     ak                    write(*,*) i
!!!$     ak                    fetxt = ' '
!!!$     ak                    errnr = 94
!!!$     ak                    goto 1000
!!!$     ak                end if

!!!$     ak                stabw = (stabw0 + stabm0/bet) * stabw
     ELSE
        IF (ldc) THEN
           IF (crtf) THEN
              READ(kanal,*,END=1001,err=1000)strnr(i),vnr(i),bet
           ELSE 
              READ(kanal,*,END=1001,err=1000)elec1,elec2,elec3,elec4,bet
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF
        ELSE
           IF (crtf) THEN
              READ(kanal,*,END=1001,err=1000)strnr(i),vnr(i),bet,pha

           ELSE
              READ(kanal,*,END=1001,err=1000)elec1,elec2,elec3,elec4,bet,pha
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF

           IF (.NOT. ldc) stabwp = ( stabpA1*bet**stabpB &
                + 1d-2*stabpA2*dabs(pha) + stabp0 ) * 1d-3
        END IF

!!!$     Ggf. Fehlermeldung
        IF (bet.LE.0d0) THEN
!!!$     ak
!!!$     ak                    write(*,*) i
           fetxt = ' '
           errnr = 94
           GOTO 1000
        END IF

        stabw = 1d-2*stabw0 + stabm0/bet

        IF ( lnse ) THEN ! add synthetic noise

           eps_r = 1d-2*nstabw0 * bet + nstabm0

           WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')rnd_r(i),eps_r,bet
! initialize
           new_pha = pha; new_bet = bet

           IF (.NOT. ldc) THEN

              eps_p = (nstabpA1*eps_r**nstabpB + 1d-2*nstabpA2*dabs(pha) &
                   + nstabp0)
              
!!$                  eps_p = (nstabpA1*bet**nstabpB + 1d-2*nstabpA2*dabs(pha) + nstabp0)
              
              WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')rnd_p(i),eps_p,pha

              new_pha = pha + rnd_p(i) * eps_p ! add noise

              WRITE(ifp2,'(G14.4)')new_pha

              IF (SIGN(pha,new_pha) /= SIGN(pha,pha)) THEN
                 
                 icount_pha = icount_pha + 1 
                 
              END IF

           END IF

           new_bet = bet + rnd_r(i) * eps_r ! add noise

           WRITE(ifp1,'(G14.4)')new_bet

           
           IF (SIGN(bet,new_bet) /= SIGN(bet,bet)) THEN
              
              icount_mag = icount_mag + 1
              
           END IF

! assign the noised values
           bet = new_bet;pha = new_pha

        END IF

     END IF

!!!$     Ggf. Fehlermeldung
     IF (bet.LE.0d0) THEN
        fetxt = ' '
        errnr = 94
        GOTO 1000
     END IF

     IF (ldc) THEN

!!!$     Phase intern auf Null setzen
        pha = 0d0
     ELSE

!!!$     Ggf. Fehlermeldung
        IF (dabs(pha).GT.1d3*pi) THEN
!!!$     ak
!!!$     ak                    write(*,*) i
           fetxt = ' '
           errnr = 95
           GOTO 1000
        END IF
     END IF

     dat(i)   = dcmplx(-dlog(bet),-pha/1d3)
     wmatdr(i) = 1d0/(stabw**2d0)
     wmatd(i) = wmatdr(i)
!!!$     ak            if (lfphai) wmatd(i)=1d0/dsqrt(stabw*stabw+stabwp*stabwp)
     IF (.NOT.ldc) THEN
        IF (lelerr) wmatd(i)=1d0/(stabw**2d0+stabwp**2d0)
        wmatdp(i)=1d0/(stabwp**2d0)
     END IF
     wdfak(i) = 1

!!!$     Stromelektroden bestimmen
     elec1 = MOD(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = MOD(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Ggf. Fehlermeldung
     IF (elec1.LT.0.OR.elec1.GT.eanz.OR. &
          elec2.LT.0.OR.elec2.GT.eanz.OR. &
          elec3.LT.0.OR.elec3.GT.eanz.OR. &
          elec4.LT.0.OR.elec4.GT.eanz) THEN
!!!$     ak
        WRITE (fetxt,'(a,I5,a)')'Electrode pair ',i,'not correct '
        errnr = 46
        GOTO 1000
     END IF
  END DO ! i=1,nanz
!!!$     ak
  IF (lindiv) THEN
     READ(kanal,*,END=1001,err=1000) stabw
     IF (stabw.LE.0d0) THEN
        fetxt = ' '
        errnr = 88
        GOTO 1000
     END IF
     DO i=1,nanz
        wmatd(i) = wmatd(i)*stabw*stabw
     END DO
  END IF
  
!!!$     'datei' schliessen
  CLOSE(kanal)
  IF ( lnse ) THEN

     IF (icount_pha > 0)THEN
        PRINT*
        PRINT*,'-- Counted ',icount_pha,&
             ' different phase signs'
!!$        fetxt = '-- check noise parameters in crt.noisemod (Phase > 0)'
!!$        errnr = 57
!!$        GOTO 1000
        PRINT*
     END IF
     IF (icount_mag > 0) THEN
        PRINT*
        PRINT*,'-- Counted ',icount_mag,&
             ' different resistance signs'
!!$        fetxt = '-- check noise parameters in crt.noisemod (Magnitude < 0)'
!!$        errnr = 57
!!$        GOTO 1000
        PRINT*
     END IF

     CLOSE (ifp1)
     IF (.NOT.ldc) CLOSE (ifp2)
     CLOSE (ifp3)
  END IF
  errnr = 0
  IF (ALLOCATED (rnd_r)) DEALLOCATE (rnd_r)
  IF (ALLOCATED (rnd_p)) DEALLOCATE (rnd_p)


  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 RETURN

1000 CLOSE(kanal)
  RETURN

1001 CLOSE(kanal)
  errnr = 2
  RETURN

END SUBROUTINE rdati
