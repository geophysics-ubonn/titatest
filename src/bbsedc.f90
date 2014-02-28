subroutine bbsedc(kanal,datei)

!!$c     Unterprogramm zur Berechnung der Summe der Sensitivitaeten aller
!!$c     Messungen (normiert)
!!$c     Berechnet nun coverage als Summe der Absolutbetraege..
!!$c     ('kpotdc' als Hilfsfeld benutzt).
!!$
!!$c     Andreas Kemna                                         02-Mar-1995
!!$c     Letzte Aenderung                   31-Mar-2010
!!$
!!$c.....................................................................

  USE alloci
  USE datmod
  USE invmod
  USE modelmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$ EIN-/AUSGABEPARAMETER:

!!!$    Kanalnummer
  INTEGER (KIND = 4) ::  kanal

!!!$    Datei
  CHARACTER (80)     ::  datei

!!!$.....................................................................

!!!$    PROGRAMMINTERNE PARAMETER:

!!!$    Hilfsvariablen
  REAL (prec) ::  xdum,ydum,dum,dum2,dummax

!!!$    Indexvariablen
  INTEGER (KIND = 4) :: i,j,k,l,i2
!!!$    Aktuelle Elementnummer
  INTEGER (KIND = 4) :: iel

!!!$    Anzahl der Knoten im aktuellen Elementtyp
  INTEGER (KIND = 4) :: nkel

!!!$.....................................................................

!!!$    Werte berechnen
  dum = 0d0

  do j = 1 , manz
     l  = j / sanz + 1
     i2 = mod(j,sanz)

     kpotdc(i2,l,1) = 0d0
     kpotdc(i2,l,2) = 0d0
     kpotdc(i2,l,3) = 0d0
     kpotdc(i2,l,4) = 0d0

     do i=1,nanz
        dum2 = SQRT(sensdc(i,j) * sensdc(i,j)) * wmatd(i) * REAL(wdfak(i))
        kpotdc(i2,l,1) = kpotdc(i2,l,1) + dum2
     end do

     dum = MAX1(dum,kpotdc(i2,l,1))
  end do
  dummax = dum
!!!$    Summe der Sensitivitaeten normieren
  if (dum.gt.1d-12) then
     do j=1,manz
        l = j/sanz + 1
        i = mod(j,sanz)
        kpotdc(i,l,1) = kpotdc(i,l,1) / dum
     end do
  end if

!!!$    Schwerpunktkoordinaten der (ggf. zusammengesetzten) Elemente bestimmen
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

        kpotdc(i2,l,2) = kpotdc(i2,l,2) + xdum/REAL(nkel)
        kpotdc(i2,l,3) = kpotdc(i2,l,3) + ydum/REAL(nkel)
        kpotdc(i2,l,4) = kpotdc(i2,l,4) + 1d0
     end do
  end do

  do j=1,manz
     l = j/sanz + 1
     i = mod(j,sanz)
     kpotdc(i,l,2) = kpotdc(i,l,2) / kpotdc(i,l,4)
     kpotdc(i,l,3) = kpotdc(i,l,3) / kpotdc(i,l,4)
  end do

!!!$    'datei' oeffnen
  fetxt = datei
  errnr = 1
  open(kanal,file=TRIM(fetxt),status='replace',err=999)
  errnr = 4

!!!$    Anzahl der Werte
  write(kanal,*,err=1000) manz

!!!$    Koordinaten und Sensitivitaetsbetraege schreiben
!!!$    (logarithmierter (Basis 10) normierter Betrag)
  do j=1,manz
     l   = j/sanz + 1
     i   = mod(j,sanz)
     dum = kpotdc(i,l,1)

     if (dum.gt.0d0) then
        write(kanal,*,err=1000) real(kpotdc(i,l,2)),&
             real(kpotdc(i,l,3)),real(LOG10(dum))
     else
        write(kanal,*,err=1000) real(kpotdc(i,l,2)),&
             real(kpotdc(i,l,3)),-1.e2
     end if
  end do

!!!$    ak
!!!$    Maximale Sensitivitaet schreiben
  write(kanal,*,err=1000)'Max:',dummax
!!!$    'datei' schliessen
  close(kanal)

  errnr = 0

  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$    Fehlermeldungen

999 return

1000 close(kanal)

END subroutine bbsedc
