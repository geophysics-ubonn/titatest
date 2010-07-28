      subroutine rdatm(kanal,datei)

c     Unterprogramm zum Einlesen der Stromwerte sowie der Elektroden-
c     kennungen aus 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   22-Feb-2006

c.....................................................................

      USE datmod
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

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4
c check whether the file format is crtomo konform or not..
      logical           ::    crtf
c
c.....................................................................

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

      ALLOCATE (strnr(nanz),strom(nanz),volt(nanz),sigmaa(nanz),
     1     kfak(nanz),vnr(nanz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation volt '
         errnr = 94
         goto 1000
      END IF
      

c     Stromelektrodennummern, Stromwerte und Spannungselektrodennummern lesen
      do i=1,nanz
         WRITE (*,'(A,I6)',ADVANCE='no')ACHAR(13)//
     1        'Getting voltage ',i
         IF (crtf) THEN
            read(kanal,*,end=1001,err=1000)
     1           strnr(i),vnr(i)
         ELSE
            read(kanal,*,end=1001,err=1000)
     1           elec1,elec2,elec3,elec4
            strnr(i) = elec1*10000 + elec2
            vnr(i)   = elec3*10000 + elec4
         END IF
c     Einheitsstrom annehmen
         strom(i) = 1d0

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
            print*,i
            fetxt = ' '
            errnr = 46
            goto 1000
         end if
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
