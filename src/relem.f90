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
  INTEGER (KIND =4)  ::    i,j,k,l,ic

!!!$     Hilfsvariable
  INTEGER (KIND =4)  ::    idum,ifln,iflnr
  LOGICAL            ::    my_check  

!!!$
  INTEGER :: ik1,ik2,jk1,jk2

!!!$ NEW rnr
  INTEGER (KIND = 4),ALLOCATABLE,DIMENSION(:) :: my_rnr

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

  my_check = .FALSE.

  do i=1,typanz
     if (typ(i).gt.10) then
        relanz = relanz + nelanz(i)
     else
        elanz  = elanz  + nelanz(i)
     end if
     my_check = my_check .OR. (typ(i) == 11)
  end do

!!!$ if all no flow boundaries, we do not have to search for a 
!!!$ average sy top...
!!$  lsytop = .NOT. my_check

!!$ get memory for the element integer field      
  ALLOCATE (nrel(elanz+relanz,smaxs),rnr(relanz),&
       my_rnr(relanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation nrel,rnr failed'
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

     end do
     idum = idum + nelanz(i)
  end do

!!!$     Zeiger auf Werte der Randelemente einlesen
  read(kanal,*,end=1001,err=1000) (rnr(i),i=1,relanz)
!!!$
!!!$ IF THIS IS NOT a pointer to the
!!!$ REGULAR ELEMENT adjacent to the border line
!!!$ THIS MAY CAUSE THE MIXED BOUNDARY TO BLOW UP
!!!$ <<< RM
  ic = 0
  DO i=1,relanz
     IF (rnr(i) > elanz.OR.rnr(i)<1) ic = ic + 1
  END DO
  IF (ic > 0) THEN
!!!$ CHECK where the BORDER ELE begin in nrel

     PRINT*,'--- Pointer of border elements are not right:'
     WRITE(*,'(a)',ADVANCE='no')'    --> rearranging'
     idum = 0;
     DO i=1,typanz
        IF (typ(i) <= 10) idum = idum + nelanz(i)
     END DO

     DO i=1,relanz
!!!$ define the node points of the border-line
        ik1 = nrel(elanz + i,1)
        ik2 = nrel(elanz + i,2)
!!!$ now search for the corresponding element
!!!$ how to do so?
!!!$ suppose we have a ordered input file, where the
!!!$ big elements come first, than we might have something like
        DO k = 1,elanz
           jk1 = 0;jk2 = 0
           DO  l = 1, selanz(1) ! TODO : fix this for multi FE
              IF (nrel(k,l) == ik1) jk1 = 1
              IF (nrel(k,l) == ik2) jk2 = 1
           END DO
           IF (jk1 == 1 .AND. jk2 == 1) THEN
              IF (lverb) PRINT*,'found border to finite element',i,k
              my_rnr(i) = k
           END IF
        END DO
     END DO
  END IF

  ic = 0
  DO i=1,relanz
     IF (my_rnr(i) > elanz.OR.my_rnr(i) < 1) ic = ic + 1
  END DO
  
  IF (ic > 0) THEN
     PRINT*,'Kann nich sein'
     fetxt = 'border to element context is wrong'
     errnr = 111
     RETURN
  ELSE
     PRINT*,'ok, copy my_rnr -> rnr'
     rnr = my_rnr
  END IF
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
