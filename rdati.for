      subroutine rdati(kanal,datei)

c     Unterprogramm zum Einlesen der Elektrodenkennungen und der Daten
c     inkl. Standardabweichungen aus 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   20-Aug-2007
      
c.....................................................................
      USE make_noise
      USE alloci, only:rnd_r,rnd_p
      USE femmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i,ifp1,ifp2

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4

c     Betrag und Phase (in mrad) der Daten
      real            * 8     bet,pha

c     Standardabweichung eines logarithmierten (!) Datums
      real            * 8     stabw
c     Error of the resistance
      real            * 8     eps_r
c     Standardabweichung der Phase
      real            * 8     stabwp,stabwb
c     Error of the phase
      real            * 8     eps_p


c     Pi
      real            * 8     pi
c    dummi
      REAL            * 8     rdum,rlev
      CHARACTER(80)     ::    csz
c check whether the file format is crtomo konform or not..
      logical           ::    crtf
c.....................................................................
      pi = dacos(-1d0)

c     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)
      errnr = 3

c     Anzahl der Messwerte lesen
      read(kanal,*,end=1001,err=1000) nanz
c check if data file format is CRTOmo konform..
      read(kanal,*,end=1001,err=1000) elec1
      BACKSPACE(kanal)

      elec3=elec1-10000 ! are we still positive?

      crtf=(elec3 > 0) ! crtomo konform?

c     Ggf. Fehlermeldung
      if (nanz.gt.nmax) then
         fetxt = ' '
         errnr = 49
         goto 1000
      end if

c     Stromelektrodennummern, Spannungselektrodennummern, Daten inkl.
c     auf 1 normierte Standardabweichungen lesen und Daten logarithmieren
      IF ( lnse ) THEN
         WRITE (*,'(A)',ADVANCE='no')ACHAR(13)//'Initializing noise'
         CALL get_unit(ifp1)
         OPEN (ifp1,FILE='inv.mynoise_rho',STATUS='replace')
         WRITE(ifp1,'(a)')'#  rnd_r'//ACHAR(9)//'eps_r'//
     1        ACHAR(9)//ACHAR(9)//'bet(old)'//ACHAR(9)//'bet(new)'
         IF (.NOT. ldc) THEN
            CALL get_unit(ifp2)
            OPEN (ifp2,FILE='inv.mynoise_pha',STATUS='replace')
            WRITE(ifp2,'(a)')'#  rnd_p'//ACHAR(9)//'eps_p'//
     1           ACHAR(9)//ACHAR(9)//'pha(old)'//ACHAR(9)//'pha(new)'
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

c$$$      WRITE (*,'(A)',ADVANCE='no')ACHAR(13)//ACHAR(9)//ACHAR(9)//
c$$$     1     ACHAR(9)//ACHAR(9)//ACHAR(9)
      do i=1,nanz
         stabwp = 0.; stabwb = 0.
         WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')ACHAR(13)//'data set ',
     1        REAL( i * (100./nanz) ),'%'
         if (lindiv) then
            if (ldc) then
               IF (crtf) THEN
                  read(kanal,*,end=1001,err=1000)
     1                 strnr(i),vnr(i),bet,stabw
               ELSE
                  read(kanal,*,end=1001,err=1000)
     1                 elec1,elec2,elec3,elec4,bet,stabw
                  strnr(i) = elec1*10000 + elec2
                  vnr(i)   = elec3*10000 + elec4
               END IF
            else
               if (lfphai) then
                  IF (crtf) THEN
                     read(kanal,*,end=1001,err=1000)
     1                    strnr(i),vnr(i),bet,pha,stabw,stabwp
                  ELSE
                     read(kanal,*,end=1001,err=1000)
     1                    elec1,elec2,elec3,elec4,bet,pha,stabw,stabwp
                     strnr(i) = elec1*10000 + elec2
                     vnr(i)   = elec3*10000 + elec4
                  END IF
c     Ggf. Fehlermeldung
                  if (stabwp.le.0d0) then
                     fetxt = ' '
                     errnr = 88
                     goto 1000
                  end if

c     ak                        stabwp = stabp0 * stabwp
                  stabwp = 1d-3*stabwp
               else
                  IF (crtf) THEN
                     read(kanal,*,end=1001,err=1000)
     1                    strnr(i),vnr(i),bet,pha,stabw
                  ELSE
                     read(kanal,*,end=1001,err=1000)
     1                    elec1,elec2,elec3,elec4,bet,pha,stabw
                     strnr(i) = elec1*10000 + elec2
                     vnr(i)   = elec3*10000 + elec4

                  END IF
               end if
            end if

c     Ggf. Fehlermeldung
            if (stabw.le.0d0) then
               fetxt = ' '
               errnr = 88
               goto 1000
            end if

c     ak                stabw = stabw0 * stabw

c     akc         Ggf. Fehlermeldung
c     ak                if (bet.le.0d0) then
c     akcak
c     ak                    write(*,*) i
c     ak                    fetxt = ' '
c     ak                    errnr = 94
c     ak                    goto 1000
c     ak                end if

c     ak                stabw = (stabw0 + stabm0/bet) * stabw
         else
            if (ldc) then
               IF (crtf) THEN
                  read(kanal,*,end=1001,err=1000)
     1                 strnr(i),vnr(i),bet
               ELSE 
                  read(kanal,*,end=1001,err=1000)
     1                 elec1,elec2,elec3,elec4,bet
                  strnr(i) = elec1*10000 + elec2
                  vnr(i)   = elec3*10000 + elec4
               END IF
            else
               IF (crtf) THEN
                  read(kanal,*,end=1001,err=1000)
     1                 strnr(i),vnr(i),bet,pha

               ELSE
                  read(kanal,*,end=1001,err=1000)
     1                 elec1,elec2,elec3,elec4,bet,pha
                  strnr(i) = elec1*10000 + elec2
                  vnr(i)   = elec3*10000 + elec4
               END IF

               if (lfphai) stabwp = ( stabpA1*bet**stabpB
     1              + 1d-2*stabpA2*dabs(pha)
     1              + stabp0 ) * 1d-3
            end if

c     Ggf. Fehlermeldung
            if (bet.le.0d0) then
c     ak
c     ak                    write(*,*) i
               fetxt = ' '
               errnr = 94
               goto 1000
            end if

            stabw = 1d-2*stabw0 + stabm0/bet

            IF ( lnse ) THEN ! add synthetic noise

               eps_r = 1d-2*nstabw0 * bet + nstabm0

               WRITE(ifp1,'(3(G14.4,1X))',ADVANCE='no')
     1              rnd_r(i),eps_r,bet

               IF (.NOT. ldc) THEN
                  
                  eps_p = (nstabpA1*eps_r**nstabpB + 
     1                 1d-2*nstabpA2*dabs(pha) + nstabp0)

c$$$                  eps_p = (nstabpA1*bet**nstabpB + 
c$$$     1                 1d-2*nstabpA2*dabs(pha) + nstabp0)

                  WRITE(ifp2,'(3(G14.4,1X))',ADVANCE='no')
     1                 rnd_p(i),eps_p,pha

                  pha = pha + rnd_p(i) * eps_p ! add noise

                  WRITE(ifp2,'(G14.4)')pha
               END IF

               bet = bet + rnd_r(i) * eps_r ! add noise

               WRITE(ifp1,'(G14.4)')bet


            END IF

         end if

c     Ggf. Fehlermeldung
         if (bet.le.0d0) then
            fetxt = ' '
            errnr = 94
            goto 1000
         end if

         if (ldc) then

c     Phase intern auf Null setzen
            pha = 0d0
         else

c     Ggf. Fehlermeldung
            if (dabs(pha).gt.1d3*pi) then
c     ak
c     ak                    write(*,*) i
               fetxt = ' '
               errnr = 95
               goto 1000
            end if
         end if

         dat(i)   = dcmplx(-dlog(bet),-pha/1d3)
         wmatdr(i) = 1d0/(stabw**2.)
         wmatd(i) = wmatdr(i)
c     ak            if (lfphai) wmatd(i)=1d0/dsqrt(stabw*stabw+stabwp*stabwp)
         IF (.NOT.ldc) THEN
            wmatd(i)=1d0/(stabw**2.+stabwp**2.)
            wmatdp(i)=1d0/(stabwp**2.)
         END IF
         wdfak(i) = 1

c     Stromelektroden bestimmen
         elec1 = mod(strnr(i),10000)
         elec2 = (strnr(i)-elec1)/10000

c     Messelektroden bestimmen
         elec3 = mod(vnr(i),10000)
         elec4 = (vnr(i)-elec3)/10000

c     Ggf. Fehlermeldung
         if (elec1.lt.0.or.elec1.gt.eanz.or.
     1        elec2.lt.0.or.elec2.gt.eanz.or.
     1        elec3.lt.0.or.elec3.gt.eanz.or.
     1        elec4.lt.0.or.elec4.gt.eanz) then
c     ak
            write(*,*) i
            fetxt = ' '
            errnr = 46
            goto 1000
         end if
      end do
c     ak
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

c     'datei' schliessen
      close(kanal)
      IF ( lnse ) THEN
         close(ifp1)
         IF (.not.ldc) close (ifp2)
      END IF
      errnr = 0
      IF (ALLOCATED (rnd_r)) DEALLOCATE (rnd_r)
      IF (ALLOCATED (rnd_p)) DEALLOCATE (rnd_p)


      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

 1001 close(kanal)
      errnr = 2
      return

      end
