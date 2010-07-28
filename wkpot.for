      subroutine wkpot(kanal,datei)

c     Unterprogramm zur Ausgabe der transformierten Potentialwerte.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   24-Jun-1997

c.....................................................................

      USE alloci
      USE electrmod
      USE elemmod
      USE wavenmod

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

c     Elektrodennummer
      integer         * 4     nelec

c     Anzahl der Knotenpunkte, an denen transformierte Potentialwerte
c     ausgegeben werden sollen
      integer         * 4     kanz

c     Nummern der Knotenpunkte, an denen transformierte Potentialwerte
c     ausgegeben werden sollen
      INTEGER(KIND = 4),DIMENSION(:),ALLOCATABLE :: knr

c     Indexvariablen
      integer         * 4     i,k

c     Hilfsvariable
      complex         * 16    dum

c.....................................................................

c     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)
      errnr = 3

c     Elektrodennummer einlesen
      read(kanal,*,end=1001,err=1000) nelec

c     Ggf. Fehlermeldung
      if (nelec.gt.eanz) then
         fetxt = ' '
         errnr = 54
         goto 1000
      end if

c     Anzahl der Knotenpunkte einlesen
      read(kanal,*,end=1001,err=1000) kanz
      
c     Ggf. Fehlermeldung
      if (kanz.gt.sanz) then
         fetxt = ' '
         errnr = 52
         goto 1000
      end if
      ALLOCATE (knr(kanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation hsens'
         errnr = 94
         GOTO 1000
      END IF
 
c     Knotennummern einlesen
      do i=1,kanz
         read(kanal,*,end=1001,err=1000) knr(i)

c     Ggf. Fehlermeldung
         if (knr(i).gt.sanz) then
            fetxt = ' '
            errnr = 53
            goto 1000
         end if
      end do

      errnr = 4

c     Entsprechenden transformierten Potentialwerte schreiben
c     (Real- und Imaginaerteil)
      do i=1,kanz
         write(kanal,*,err=1000)
         write(kanal,*,err=1000) knr(i)

         do k=1,kwnanz
            dum = kpot(knr(i),nelec,k)
            write(kanal,*,err=1000) real(dble(dum)),
     1           real(dimag(dum))
         end do
      end do

c     'datei' schliessen
      close(kanal)

      errnr = 0
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
