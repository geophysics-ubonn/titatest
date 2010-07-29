      subroutine welectr(kanal,datei)

c     Unterprogramm zum Schreiben der Elektrodenverteilung in 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   21-Jan-2003

c.....................................................................

      USE electrmod
      
      IMPLICIT none

      INCLUDE 'err.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i

c.....................................................................

c     'datei' oeffnen
      fetxt = datei

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=1000)

      errnr = 4

c     Anzahl der Elektroden schreiben
      write(kanal,*,err=1000) eanz

c     Knotennummern der Elektroden schreiben
      do i=1,eanz
         write(kanal,*,err=1000) enr(i)
      end do

c     'datei' schliessen
      close(kanal)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 1000 return

      end
