subroutine bpot(kanal,datei)

!!!$     Unterprogramm zur Berechnung der Potentialwerte
!!!$     (beachte: Potentialwerte wurden fuer Einheitsstrom berechnet).

!!!$     Andreas Kemna                                            17-Aug-1994
!!!$     Letzte Aenderung   24-Jun-1997

!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4)  ::     kanal

!!!$     Datei
  CHARACTER (80)      ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j

!!!$     Elektrodennummern
  INTEGER (KIND = 4)  ::     elec1,elec2

!!!$.....................................................................

  do i=1,nanz

!!!$     Stromelektroden bestimmen
     elec1 = mod(strnr(i),10000)
     elec2 = (strnr(i)-elec1)/10000

     do j=1,sanz

!!!$     (beachte: Faktoren '.../2d0' (-> Potentialwerte fuer Einheitsstrom)
!!!$     und '...*2d0' (-> Ruecktransformation) kuerzen sich weg !)
        if (elec1.eq.0) then
           pot(j) = hpot(j,elec2) * dcmplx(strom(i))
        else if (elec2.eq.0) then
           pot(j) = -hpot(j,elec1) * dcmplx(strom(i))
        else
           pot(j) = (hpot(j,elec2)-hpot(j,elec1)) * dcmplx(strom(i))
        end if
     end do

!!!$     Potentialwerte ausgeben
     call wpot(kanal,datei,i)
     if (errnr.ne.0) goto 1000
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine bpot
