      subroutine relectr(kanal,datei)
      
!!!$     Unterprogramm zum Einlesen der Elektrodenverteilung aus 'datei'.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   24-Oct-1996

!!!$.....................................................................

      USE electrmod
      USE elemmod
      USE errmod

      IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
      integer         * 4     kanal

!!!$     Datei
      character       * 80    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
      integer         * 4     i,ifp

!!!$.....................................................................

!!!$     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)
      CALL get_unit(ifp)

      OPEN (ifp,FILE='inv.elecpositions',STATUS='replace')

      errnr = 3

!!!$     Anzahl der Elektroden einlesen
      read(kanal,*,end=1001,err=1000) eanz
!!!$ memory allocation
      ALLOCATE (enr(eanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation enr'
         errnr = 94
         goto 1000
      END IF

      WRITE (ifp,*)eanz

!!!$     Knotennummern der Elektroden einlesen
      do i=1,eanz
         read(kanal,*,end=1001,err=1000) enr(i)
         WRITE (ifp,*)sx(snr(enr(i))),sy(snr(enr(i)))
!!!$     Ggf. Fehlermeldung
         if (enr(i).gt.sanz) then
            fetxt = ' '
            errnr = 29
            goto 1000
         end if
      end do

!!!$     'datei' schliessen
      CLOSE (ifp)
      close(kanal)

      errnr = 0
      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

 1001 close(kanal)
      errnr = 2
      return

      end
