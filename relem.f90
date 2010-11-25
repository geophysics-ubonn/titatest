subroutine relem(kanal,datei)

!!!$     Unterprogramm zum Einlesen der FEM-Parameter aus 'datei'.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   24-Oct-1996

!!!$.....................................................................

  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT none


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4) ::    kanal

!!!$     Datei
  CHARACTER (80)     ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariablen
  INTEGER (KIND =4)  ::    i,j,k

!!!$     Hilfsvariable
  INTEGER (KIND =4)  ::    idum,ifln,iflnr

!!!$.....................................................................

!!!$     'datei' oeffnen
  fetxt = datei

  errnr = 1
  open(kanal,file=TRIM(fetxt),status='old',err=999)

  errnr = 3

!!!$     Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
!!!$     sowie Bandbreite der Gesamtsteifigkeitsmatrix einlesen
  read(kanal,*,end=1001,err=1000) sanz,typanz,mb

!!$ now get some memory for the fields..
!!$ first the sanz fields
  ALLOCATE (sx(sanz),sy(sanz),snr(sanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation sx failed'
     errnr = 94
     GOTO 999
  END IF

  ALLOCATE (typ(typanz),nelanz(typanz),selanz(typanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation selanz failed'
     errnr = 94
     GOTO 999
  END IF

!!!$     Elementtypen, Anzahl der Elemente eines bestimmten Typs sowie
!!!$     Anzahl der Knoten in einem Elementtyp einlesen
  read(kanal,*,end=1001,err=1000)(typ(i),nelanz(i),selanz(i),i=1,typanz)

!!$ set number of node points for regular elements
  smaxs = MAXVAL(selanz)

!!!$     Anzahl der Elemente (ohne Randelemente) und Anzahl der Randelemente
!!!$     bestimmen
  relanz = 0
  elanz  = 0

  do i=1,typanz
     if (typ(i).gt.10) then
        relanz = relanz + nelanz(i)
     else
        elanz  = elanz  + nelanz(i)
     end if
  end do

!!$ get memory for the element integer field      
  ALLOCATE (nrel(elanz+relanz,smaxs),rnr(relanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation nrel failed'
     errnr = 94
     GOTO 999
  END IF
!!$ get memory for the regular element midpoint coordinates
  ALLOCATE (espx(elanz),espy(elanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation espx failed'
     errnr = 94
     GOTO 999
  END IF
  espx = 0.;espy = 0.
!!!$     Zeiger auf Koordinaten, x-Koordinaten sowie y-Koordinaten der Knoten
!!!$     einlesen
  read(kanal,*,end=1001,err=1000) (snr(i),sx(i),sy(i),i=1,sanz)
!!!$     Knotennummern der Elemente einlesen
  idum = 0;ifln = 0;iflnr = 0
  do i=1,typanz
     do j=1,nelanz(i)
        read(kanal,*,end=1001,err=1000)(nrel(idum+j,k),k=1,selanz(i))

        IF (typ(i) < 10) THEN ! set midpoints

           ifln = ifln + 1

           DO k = 1,selanz(i)
              espx(ifln) = espx(ifln) + sx(snr(nrel(idum+j,k)))
              espy(ifln) = espy(ifln) + sy(snr(nrel(idum+j,k)))
           END DO

           espx(ifln) = espx(ifln) / selanz(i)
           espy(ifln) = espy(ifln) / selanz(i)

        END IF

!!!$            IF (typ(i) > 10) THEN ! randele zeiger kann man auch so belegen
!!!$               iflnr = iflnr + 1
!!!$               IF (iflnr > relanz) THEN
!!!$                  fetxt = 'relem:: iflnr > relanz!'
!!!$                  errnr = 9
!!!$                  goto 1000
!!!$               END IF                  
!!!$               rnr(iflnr) = nrel(idum+j,1)
!!!$            END IF
     end do
     idum = idum + nelanz(i)
  end do

!!!$     Zeiger auf Werte der Randelemente einlesen
  read(kanal,*,end=1001,err=1000) (rnr(i),i=1,relanz)

!!!$     'datei' schliessen
  close(kanal)

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 return

1000 close(kanal)
  return

1001 close(kanal)
  errnr = 2
  return

end subroutine relem
