subroutine bbsens(kanal,datei)

!!!$     Unterprogramm zur Berechnung der Summe der Sensitivitaeten 
!!!$     aller Messungen (normiert)
!!!$     Berechnet nun coverage als Summe der Absolutbetraege..
!!!$     ('kpot' als Hilfsfeld benutzt).

!!!$     Andreas Kemna                                   02-Mar-1995
!!!$     Letzte Aenderung              31-Mar-2010

!!!$....................................................................

  USE alloci
  USE datmod
  USE invmod
  USE modelmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4)  ::     kanal

!!!$     Datei
  CHARACTER (80)      ::  datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  REAL (KIND(0D0))    ::     xdum,ydum,dum,dum2,dummax

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j,k,l,i2

!!!$     Aktuelle Elementnummer
  INTEGER (KIND = 4)  ::     iel

!!!$     Anzahl der Knoten im aktuellen Elementtyp
  INTEGER (KIND = 4)  ::     nkel

!!!$.....................................................................

!!!$     Werte berechnen
  dum = 0d0

  do j=1,manz
     l  = j/sanz + 1
     i2 = mod(j,sanz)

     kpot(i2,l,1) = dcmplx(0d0)
     kpot(i2,l,2) = dcmplx(0d0)
     kpot(i2,l,3) = dcmplx(0d0)
     kpot(i2,l,4) = dcmplx(0d0)

     do i=1,nanz            ! coverage -> L1 Norm
        dum2 = SQRT(dconjg(sens(i,j)) * sens(i,j) * &
             dcmplx(wmatd(i)*dble(wdfak(i))))
        kpot(i2,l,1) = kpot(i2,l,1) + dum2
     end do

     dum = dmax1(dum,dble(kpot(i2,l,1)))
  end do
  dummax = dum
!!!$     Summe der Sensitivitaeten normieren
  if (dum.gt.1d-12) then
     do j=1,manz
        l = j/sanz + 1
        i = mod(j,sanz)
        kpot(i,l,1) = kpot(i,l,1) / dcmplx(dum)
     end do
  end if

!!!$     Schwerpunktkoordinaten der (ggf. zusammengesetzten) Elemente bestimmen
  iel = 0

  do i=1,typanz
     if (typ(i).gt.10) EXIT

     nkel = selanz(i)

     do j=1,nelanz(i)
        iel  = iel + 1

        xdum = 0d0
        ydum = 0d0

        do k=1,nkel
           xdum = xdum + sx(snr(nrel(iel,k)))
           ydum = ydum + sy(snr(nrel(iel,k)))
        end do

        l  = mnr(iel)/sanz + 1
        i2 = mod(mnr(iel),sanz)

        kpot(i2,l,2) = kpot(i2,l,2) + dcmplx(xdum/dble(nkel))
        kpot(i2,l,3) = kpot(i2,l,3) + dcmplx(ydum/dble(nkel))
        kpot(i2,l,4) = kpot(i2,l,4) + dcmplx(1d0)
     end do
  end do

  do j=1,manz
     l = j/sanz + 1
     i = mod(j,sanz)
     kpot(i,l,2) = kpot(i,l,2) / kpot(i,l,4)
     kpot(i,l,3) = kpot(i,l,3) / kpot(i,l,4)
  end do

!!!$     'datei' oeffnen
  fetxt = datei
  errnr = 1
  open(kanal,file=fetxt,status='replace',err=999)
  errnr = 4

!!!$     ak
!!!$     Maximale Sensitivitaet schreiben
  write(kanal,*,err=1000) manz

!!!$     Koordinaten und Sensitivitaetsbetraege schreiben
!!!$     (logarithmierter (Basis 10) normierter Betrag)
  do j=1,manz
     l   = j/sanz + 1
     i   = mod(j,sanz)
     dum = dble(kpot(i,l,1))

     if (dum.gt.0d0) then
        write(kanal,*,err=1000) real(dble(kpot(i,l,2))),&
             real(dble(kpot(i,l,3))),real(dlog10(dum))
     else
        write(kanal,*,err=1000) real(dble(kpot(i,l,2))),&
             real(dble(kpot(i,l,3))),-1.e2
     end if
  end do

  write(kanal,*,err=1000)'Max:',dummax
!!!$     'datei' schliessen
  close(kanal)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

end subroutine bbsens
