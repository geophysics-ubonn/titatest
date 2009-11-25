      subroutine wdatm(kanal,datei)

c     Unterprogramm zum Schreiben der Strom- und Spannungs- bzw.
c     scheinbaren Widerstandswerte sowie der Elektrodenkennungen
c     in 'datei'.

c     Andreas Kemna                                            22-Oct-1993
c     Letzte Aenderung   10-Mar-2007

c.....................................................................

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
      INCLUDE 'dat.fin'

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

c     Pi
      real            * 8     pi

c     Hilfsvariablen
      real            * 8     bet,pha,npi
      integer         * 4     ie1,ie2
c     Standartabweichung.-..
      real            * 4     stab
c.....................................................................

 100  FORMAT(2(4X,I6),2(1x,G14.7))
      pi = dacos(-1d0)

c     'datei' oeffnen
      fetxt = datei
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

c     Anzahl der Messwerte schreiben
      write(kanal,*,err=1000) nanz

c     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
c     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
      stab=5.0
      do i=1,nanz
         bet = cdabs(sigmaa(i))
         pha = datan2(dimag(sigmaa(i)),dble(sigmaa(i)))

c     ak
c     Ggf. Polaritaet vertauschen
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
c     1                  strnr(i),vnr(i),real(bet),real(1d3*pha),
c     1                  real(kfak(i))
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

      end
