      subroutine wdatm(kanal,datei)

!!!$     Unterprogramm zum Schreiben der Strom- und Spannungs- bzw.
!!!$     scheinbaren Widerstandswerte sowie der Elektrodenkennungen
!!!$     in 'datei'.

!!!$     Andreas Kemna                                            22-Oct-1993
!!!$     Letzte Aenderung   10-Mar-2007

!!!$.....................................................................

      USE datmod
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
      integer         * 4     i

!!!$     Pi
      real            * 8     pi

!!!$     Hilfsvariablen
      real            * 8     bet,pha,npi
      integer         * 4     ie1,ie2
!!!$     Standartabweichung.-..
      real            * 4     stab
!!!$.....................................................................

 100  FORMAT(2(4X,I6),2(1x,G14.7))
      pi = dacos(-1d0)

!!!$     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

!!!$     Anzahl der Messwerte schreiben
      write(kanal,*,err=1000) nanz

!!!$     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
!!!$     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
      stab=5.0
      do i=1,nanz
         bet = cdabs(sigmaa(i))
         pha = datan2(dimag(sigmaa(i)),dble(sigmaa(i)))

!!!$     ak
!!!$     Ggf. Polaritaet vertauschen
         npi = dnint(pha/pi)*pi
         if (dabs(npi).gt.1d-12) then
            pha    = pha-npi
            ie1    = mod(vnr(i),10000)
            ie2    = (vnr(i)-ie1)/10000
            vnr(i) = ie1*10000+ie2
         end if

c$$$  write(kanal,*,err=1000)
c$$$  1                  strnr(i),vnr(i),real(bet),real(1d3*pha)
         write(kanal,100,err=1000)
     1        strnr(i),vnr(i),real(bet),real(1d3*pha)
!!!$     1                  strnr(i),vnr(i),real(bet),real(1d3*pha),
!!!$     1                  real(kfak(i))
      end do

!!!$     'datei' schliessen
      close(kanal)

      errnr = 0
      return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end
