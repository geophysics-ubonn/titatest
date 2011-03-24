subroutine rdati(kanal,datei)

!!!$     Unterprogramm zum Einlesen der Elektrodenkennungen und der Daten
!!!$     inkl. Standardabweichungen aus 'datei'.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   20-Aug-2007

!!!$.....................................................................
  USE make_noise
  USE alloci, only:rnd_r,rnd_p
  USE femmod
  USE datmod
  USE invmod
  USE electrmod
  USE errmod
  USE konvmod

  IMPLICIT none


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
  open(kanal,file=TRIM(fetxt),status='old',err=999)
  errnr = 3

!!!$     Anzahl der Messwerte lesen
  read(kanal,*,end=1001,err=1000) nanz
!!!$c check if data file format is CRTOmo konform..
  read(kanal,*,end=1001,err=1000) elec1
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
     goto 1000
  END IF

!!!$     Stromelektrodennummern, Spannungselektrodennummern, Daten inkl.
!!!$     auf 1 normierte Standardabweichungen lesen und Daten logarithmieren
  IF ( lnse ) THEN
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


  do i=1,nanz
     stabwp = 0.; stabwb = 0.
     IF (lverb) WRITE (*,'(A,t70,F6.2,A)',ADVANCE='no')ACHAR(13)//&
          'data set ',REAL( i * (100./nanz) ),'%'
     if (lindiv) then
        if (ldc) then
           IF (crtf) THEN
              read(kanal,*,end=1001,err=1000)strnr(i),vnr(i),bet,stabw
           ELSE
              read(kanal,*,end=1001,err=1000)elec1,elec2,elec3,elec4,bet,&
                   stabw
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF
        else
           if (lfphai) then
              IF (crtf) THEN
                 read(kanal,*,end=1001,err=1000)strnr(i),vnr(i),bet,pha,&
                      stabw,stabwp
              ELSE
                 read(kanal,*,end=1001,err=1000)elec1,elec2,elec3,elec4,&
                      bet,pha,stabw,stabwp
                 strnr(i) = elec1*10000 + elec2
                 vnr(i)   = elec3*10000 + elec4
              END IF
!!!$     Ggf. Fehlermeldung
              if (stabwp.le.0d0) then
                 fetxt = ' '
                 errnr = 88
                 goto 1000
              end if

!!!$     ak                        stabwp = stabp0 * stabwp
              stabwp = 1d-3*stabwp
           else
              IF (crtf) THEN
                 read(kanal,*,end=1001,err=1000)strnr(i),vnr(i),bet,pha,&
                      stabw
              ELSE
                 read(kanal,*,end=1001,err=1000)elec1,elec2,elec3,elec4,&
                      bet,pha,stabw
                 strnr(i) = elec1*10000 + elec2
                 vnr(i)   = elec3*10000 + elec4

              END IF
           end if
        end if

!!!$     Ggf. Fehlermeldung
        if (stabw.le.0d0) then
           fetxt = ' '
           errnr = 88
           goto 1000
        end if

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
     else
        if (ldc) then
           IF (crtf) THEN
              read(kanal,*,end=1001,err=1000)strnr(i),vnr(i),bet
           ELSE 
              read(kanal,*,end=1001,err=1000)elec1,elec2,elec3,elec4,bet
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF
        else
           IF (crtf) THEN
              read(kanal,*,end=1001,err=1000)strnr(i),vnr(i),bet,pha

           ELSE
              read(kanal,*,end=1001,err=1000)elec1,elec2,elec3,elec4,bet,pha
              strnr(i) = elec1*10000 + elec2
              vnr(i)   = elec3*10000 + elec4
           END IF

           if (.NOT. ldc) stabwp = ( stabpA1*bet**stabpB &
                + 1d-2*stabpA2*dabs(pha) + stabp0 ) * 1d-3
        end if

!!!$     Ggf. Fehlermeldung
        if (bet.le.0d0) then
!!!$     ak
!!!$     ak                    write(*,*) i
           fetxt = ' '
           errnr = 94
           goto 1000
        end if

        stabw = 1d-2*stabw0 + stabm0/bet

        IF ( lnse ) THEN ! add synthetic noise

           eps_r = 1d-2*nstabw0 * bet + nstabm0

           WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')rnd_r(i),eps_r,bet

           IF (.NOT. ldc) THEN

              eps_p = (nstabpA1*eps_r**nstabpB + 1d-2*nstabpA2*dabs(pha) &
                   + nstabp0)
              
!!$                  eps_p = (nstabpA1*bet**nstabpB + 1d-2*nstabpA2*dabs(pha) + nstabp0)
              
              WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')rnd_p(i),eps_p,pha

              pha = pha + rnd_p(i) * eps_p ! add noise

              WRITE(ifp2,'(G14.4)')pha
           END IF

           bet = bet + rnd_r(i) * eps_r ! add noise

           WRITE(ifp1,'(G14.4)')bet
           ! write out full noisy data as measured voltages..
           WRITE (ifp3,*) strnr(i),vnr(i),bet,pha
        END IF

     end if

!!!$     Ggf. Fehlermeldung
     if (bet.le.0d0) then
        fetxt = ' '
        errnr = 94
        goto 1000
     end if

     if (ldc) then

!!!$     Phase intern auf Null setzen
        pha = 0d0
     else

!!!$     Ggf. Fehlermeldung
        if (dabs(pha).gt.1d3*pi) then
!!!$     ak
!!!$     ak                    write(*,*) i
           fetxt = ' '
           errnr = 95
           goto 1000
        end if
     end if

     dat(i)   = dcmplx(-dlog(bet),-pha/1d3)
     wmatdr(i) = 1d0/(stabw**2.)
     wmatd(i) = wmatdr(i)
!!!$     ak            if (lfphai) wmatd(i)=1d0/dsqrt(stabw*stabw+stabwp*stabwp)
     IF (.NOT.ldc) THEN
        IF (lelerr) wmatd(i)=1d0/(stabw**2.+stabwp**2.)
        wmatdp(i)=1d0/(stabwp**2.)
     END IF
     wdfak(i) = 1

!!!$     Stromelektroden bestimmen
     elec1 = mod(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

!!!$     Messelektroden bestimmen
     elec3 = mod(vnr(i),10000)
     elec4 = (vnr(i)-elec3)/10000

!!!$     Ggf. Fehlermeldung
     if (elec1.lt.0.or.elec1.gt.eanz.or. &
          elec2.lt.0.or.elec2.gt.eanz.or. &
          elec3.lt.0.or.elec3.gt.eanz.or. &
          elec4.lt.0.or.elec4.gt.eanz) then
!!!$     ak
        WRITE (fetxt,'(a,I5,a)')'Electrode pair ',i,'not correct '
        errnr = 46
        goto 1000
     end if
  end do
!!!$     ak
  if (lindiv) then
     read(kanal,*,end=1001,err=1000) stabw
     if (stabw.le.0d0) then
        fetxt = ' '
        errnr = 88
        goto 1000
     end if
     do i=1,nanz
        wmatd(i) = wmatd(i)*stabw*stabw
     end do
  end if

!!!$     'datei' schliessen
  close(kanal)
  IF ( lnse ) THEN
     CLOSE (ifp1)
     IF (.not.ldc) CLOSE (ifp2)
     CLOSE (ifp3)
  END IF
  errnr = 0
  IF (ALLOCATED (rnd_r)) DEALLOCATE (rnd_r)
  IF (ALLOCATED (rnd_p)) DEALLOCATE (rnd_p)


  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

1001 close(kanal)
  errnr = 2
  return

end subroutine rdati
