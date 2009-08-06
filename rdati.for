      subroutine rdati(kanal,datei)

c     Unterprogramm zum Einlesen der Elektrodenkennungen und der Daten
c     inkl. Standardabweichungen aus 'datei'.

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   20-Aug-2007
      
c.....................................................................
      USE make_noise

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'electr.fin'
      INCLUDE 'fem.fin'
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
      integer         * 4     i

c     Elektrodennummern
      integer         * 4     elec1,elec2,
     1     elec3,elec4

c     Betrag und Phase (in mrad) der Daten
      real            * 8     bet,pha

c     Standardabweichung eines logarithmierten (!) Datums
      real            * 8     stabw

c     Standardabweichung der Phase
      real            * 8     stabwp

c     Pi
      real            * 8     pi

c.....................................................................

      pi = dacos(-1d0)

c     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='old',err=999)
      errnr = 3

c     Anzahl der Messwerte lesen
      read(kanal,*,end=1001,err=1000) nanz

c     Ggf. Fehlermeldung
      if (nanz.gt.nmax) then
         fetxt = ' '
         errnr = 49
         goto 1000
      end if

c     Stromelektrodennummern, Spannungselektrodennummern, Daten inkl.
c     auf 1 normierte Standardabweichungen lesen und Daten logarithmieren
      do i=1,nanz
         WRITE (*,'(A,I9)',advance='no')'data set '//ACHAR(13),i
         if (lindiv) then
            if (ldc) then
               read(kanal,*,end=1001,err=1000)
c     ro,ERT2003
     1              strnr(i),vnr(i),bet,stabw
c     ro,ERT2003     1                         strnr(i),vnr(i),bet,pha,stabw
            else
               if (lfphai) then
                  read(kanal,*,end=1001,err=1000)
     1                 strnr(i),vnr(i),bet,pha,stabw,stabwp

c     Ggf. Fehlermeldung
                  if (stabwp.le.0d0) then
                     fetxt = ' '
                     errnr = 88
                     goto 1000
                  end if

c     ak                        stabwp = stabp0 * stabwp
                  stabwp = 1d-3*stabwp
               else
                  read(kanal,*,end=1001,err=1000)
     1                 strnr(i),vnr(i),bet,pha,stabw
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
               read(kanal,*,end=1001,err=1000)
     1              strnr(i),vnr(i),bet
c     ak Inga
c     ak     1                         elec1,elec2,elec3,elec4,bet
c     ak                    strnr(i) = elec1*10000 + elec2
c     ak                    vnr(i)   = elec3*10000 + elec4

            else
               read(kanal,*,end=1001,err=1000)
     1              strnr(i),vnr(i),bet,pha

               if (lfphai) stabwp = 1d-3*stabpA1*bet**stabpB
     1              + 1d-2*stabpA2*1d-3*dabs(pha)
     1              + 1d-3*stabp0
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
         wmatd(i) = 1d0/(stabw*stabw)
c     ak            if (lfphai) wmatd(i)=1d0/dsqrt(stabw*stabw+stabwp*stabwp)
         if (lfphai) wmatdp(i)=1d0/(stabwp*stabwp)
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
