subroutine wdatm(kanal,datei)

!!!$     Unterprogramm zum Schreiben der Strom- und Spannungs- bzw.
!!!$     scheinbaren Widerstandswerte sowie der Elektrodenkennungen
!!!$     in 'datei'.

!!!$     Andreas Kemna                                            22-Oct-1993
!!!$     Letzte Aenderung   10-Mar-2007

!!!$.....................................................................
use alloci, only:prec
  USE datmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND=4) ::   kanal

!!!$     Datei
  CHARACTER (80)   ::  datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariable
  INTEGER (KIND=4) ::    i

!!!$     Pi
  REAL(prec)  ::    pi

!!!$     Hilfsvariablen
  REAL(prec)  ::   bet,pha,npi
  INTEGER (KIND=4) ::    ie1,ie2
!!!$     Standartabweichung.-..
  REAL(prec)  ::    stab
!!!$.....................................................................

  pi = dacos(-1d0)

!!!$     'datei' oeffnen
  fetxt = datei
  errnr = 1
  open(kanal,file=TRIM(fetxt),status='replace',err=999)
  errnr = 4

!!!$     Anzahl der Messwerte schreiben
  write(kanal,*,err=1000) nanz

!!!$     Stromelektrodennummern, Spannungselektrodennummern und scheinbare
!!!$     Widerstandswerte (Betrag und Phase (in mrad)) schreiben
  stab=5.0
  do i=1,nanz
     bet = ABS(sigmaa(i))
     pha = ATAN2(aimag(sigmaa(i)),REAL(sigmaa(i)))

!!!$     ak
!!!$     Ggf. Polaritaet vertauschen
     npi = nint(pha/pi)*pi
     if (ABS(npi).gt.1d-12) then
        pha    = pha-npi
        ie1    = mod(vnr(i),10000)
        ie2    = (vnr(i)-ie1)/10000
        vnr(i) = ie1*10000+ie2
     end if

!!!$  write(kanal,*,err=1000)
!!!$  1                  strnr(i),vnr(i),real(bet),real(1d3*pha)
     write(kanal,*,err=1000)strnr(i),vnr(i),real(bet),real(1d3*pha)
!!!$     1                  strnr(i),vnr(i),real(bet),real(1d3*pha),
!!!$     1                  real(kfak(i))
  end do
write(*,'(a,I7,a)') ' wrote',nanz,' voltages'
!!!$     'datei' schliessen
  close(kanal)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

end subroutine wdatm
