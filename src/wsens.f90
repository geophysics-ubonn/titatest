subroutine wsens(kanal,datei)

!!!$     Unterprogramm zum Schreiben der Sensitivitaeten aller Messungen.

!!!$     Andreas Kemna                                            20-Apr-1995
!!!$     Letzte Aenderung   10-Mar-2007

!!!$.....................................................................

  USE datmod
  USE alloci
  USE modelmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND=4) ::       kanal

!!!$     Datei
  CHARACTER (80)   ::     datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsfunctions
  CHARACTER (12)   ::    intcha
  CHARACTER (80)   ::    filpat

!!!$     Hilfsvariablen
  INTEGER (KIND=4) ::     idum,idum2,lnanz
  CHARACTER (80)   ::    htxt
  CHARACTER (12)   ::    htxt2

  COMPLEX(prec) ::    summe
!!!$     Indexvariablen
  INTEGER (KIND=4) ::     i,j

!!!$     (Back-) Slash
  CHARACTER (1)   ::     slash

!!!$.....................................................................

!!!$     Slash
  slash = '/'

!!!$     'datei' modifizieren
  lnanz = int(log10(real(nanz)))+1
  htxt  = filpat(datei,idum2,1,slash(1:1))
  idum  = idum2+index(datei(idum2+1:80),'.')-1

  if ((idum-idum2-1).gt.(8-lnanz)) then
     fetxt = datei
     errnr = 15
     goto 999
  end if

!!!$     Ausgabeschleife
  do i=1,nanz

!!!$     'datei' bestimmen
     htxt2 = intcha(i,lnanz)
     htxt  = datei(1:idum)//htxt2(1:lnanz)//datei(idum+1:idum+4)

!!!$     'datei' oeffnen
     fetxt = htxt
     errnr = 1
     open(kanal,file=TRIM(fetxt),status='replace',err=999)
     errnr = 4

!!!$     Summe der Sensitivitaeten berechnen und ausgeben
!!!$     (Betrag und Phase (in mrad))
!!!$     ak            sensmax = 0d0
     summe   = CMPLX(0d0)

     do j=1,manz
!!!$     ak                sensmax = MAX1(sensmax,ABS(sens(i,j)))
        summe   = summe + sens(i,j)
     end do

     if (ABS(summe).gt.0d0) then
        write(kanal,*,err=1000) real(ABS(summe)),&
             real(1d3*ATAN2(aimag(summe),real(summe)))
     else
        write(kanal,*,err=1000) 0.,0.
     end if

!!!$     Koordinaten und Sensitivitaeten schreiben
!!!$     (Real- und aimaginaerteil)
     do j=1,manz
        write(kanal,*,err=1000)real(espx(j)),real(espy(j)),&
             real(REAL(sens(i,j))),real(aimag(sens(i,j)))
     end do

!!!$     'datei' schliessen
     close(kanal)
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

end subroutine wsens
