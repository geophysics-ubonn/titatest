subroutine elem8(kelmas,kelve,kwert,smaxs)

!!!$     Unterprogramm liefert die Elementmatrix 'kelmas(4,4)' und den Element-
!!!$     vektor 'kelve(4)' fuer ein zusammengesetztes Viereckelement
!!!$     (vier Teildreiecke mit linearem Ansatz) ( Elementtyp Nr.8 ).

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   26-Jan-1998

!!!$.....................................................................
use alloci, only: prec
  USE elemmod,only:xk,yk
  USE errmod

  IMPLICIT none


!!!$.....................................................................
  INTEGER (KIND = 4)                     :: smaxs
!!!$     Elementmatrix nach der Kondensation
  REAL(prec),DIMENSION(smaxs,smaxs) :: kelmas

!!!$     Elementvektor nach der Kondensation
  REAL(prec),DIMENSION(smaxs)       :: kelve

!!!$     Wellenzahlwert
  REAL(prec)                        :: kwert

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfskoordinaten
  REAL (prec)    ::     xkdum(5),ykdum(5)

!!!$     Beteiligte Knoten des zusammengesetzten Elements am Teildreieck
  INTEGER (KIND = 4)  ::     ik(3)

!!!$     Elementmatrizen vor der Kondensation
  REAL (prec)    ::     elmas(5,5),elmam(5,5)

!!!$     Elementvektor vor der Kondensation
  REAL (prec)    ::     elve(5)

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j,k

!!!$     Grundelementmatrizen
  INTEGER (KIND = 4)  ::     s1(3,3),s2(3,3),s3(3,3),s4(3,3)

!!!$     Grundelementvektor
  INTEGER (KIND = 4)  ::     sb(3)

!!!$     Hilfsvariablen
  REAL (prec)    ::     x21,x31,y21,y31,det,a,b,c,a1,a2

!!!$.....................................................................

  data s1/1,-1,0,-1,1,0,0,0,0/
  data s2/2,-1,-1,-1,0,1,-1,1,0/
  data s3/1,0,-1,0,0,0,-1,0,1/
  data s4/2,1,1,1,2,1,1,1,2/
  data sb/1,1,1/

!!!$     Elementmatrix bzw. -vektor des zusammengesetzten Elements gleich
!!!$     Null setzen
  do i=1,5
     elve(i) = 0d0

     do j=1,5
        elmas(i,j) = 0d0
        elmam(i,j) = 0d0
     end do
  end do

!!!$     Hilfskoordinaten bestimmen
  a1 = xk(3)-xk(1)
  a2 = xk(4)-xk(2)

!!!$     Ggf. Fehlermeldung
  if (ABS(a1).le.1d-12.or.ABS(a2).le.1d-12) then
     fetxt = ' '
!!!$     ak            write(fetxt(1:20),'(g20.5)') ABS(a1)
!!!$     ak            write(fetxt(26:45),'(g20.5)') ABS(a2)
     errnr = 26
     goto 1000
  end if

  a1 = (yk(3)-yk(1))/a1
  a2 = (yk(4)-yk(2))/a2

!!!$     Ggf. Fehlermeldung
  if (ABS(a1-a2).le.1d-12) then
     fetxt = ' '
!!!$     ak            write(fetxt(1:20),'(g20.5)') ABS(a1-a2)
     errnr = 26
     goto 1000
  end if

  xkdum(5) = (a1*xk(1)-yk(1)-a2*xk(2)+yk(2))/(a1-a2)
  ykdum(5) = yk(1)+a1*(xkdum(5)-xk(1))

!!!$     Restlichen Koordinaten umspeichern
  do i=1,4
     xkdum(i) = xk(i)
     ykdum(i) = yk(i)
  end do

!!!$     Beitraege der Teilelemente aufaddieren
  do k=1,4

!!!$     Beteiligte Knotenpunkte definieren
     if (k.eq.1) then
        ik(1) = 1
        ik(2) = 2
        ik(3) = 5
     else if (k.eq.2) then
        ik(1) = 2
        ik(2) = 3
        ik(3) = 5
     else if (k.eq.3) then
        ik(1) = 5
        ik(2) = 3
        ik(3) = 4
     else if (k.eq.4) then
        ik(1) = 1
        ik(2) = 5
        ik(3) = 4
     end if

     x21 = xkdum(ik(2)) - xkdum(ik(1))
     x31 = xkdum(ik(3)) - xkdum(ik(1))
     y21 = ykdum(ik(2)) - ykdum(ik(1))
     y31 = ykdum(ik(3)) - ykdum(ik(1))
     det = x21*y31 - x31*y21

!!!$     Ggf. Fehlermeldung
     if (det.le.1d-12) then
        fetxt = ' '
        write(fetxt(1:20),'(g20.5)') det
        errnr = 26
        goto 1000
     end if

     a =   (x31*x31 + y31*y31) / det
     b = - (x31*x21 + y31*y21) / det
     c =   (x21*x21 + y21*y21) / det

     do i=1,3
        elve(ik(i)) = elve(ik(i)) + det * REAL(sb(i)) / 6d0

        do j=1,3
           elmas(ik(i),ik(j)) = elmas(ik(i),ik(j)) + &
                (a*REAL(s1(i,j)) + b*REAL(s2(i,j)) + &
                c*REAL(s3(i,j))) / 2d0
           elmam(ik(i),ik(j)) = elmam(ik(i),ik(j)) + &
                det * REAL(s4(i,j)) / 2.4d1
        end do
     end do
  end do

!!!$     Kondensationsschritt
  do i=1,5
     do j=1,5
        elmas(i,j) = elmas(i,j) + kwert*kwert*elmam(i,j)
     end do
  end do

  do i=1,4
     kelve(i) = elve(i) - elmas(i,5)*elve(5)/elmas(5,5)

     do j=1,4
        kelmas(i,j) = elmas(i,j) - elmas(i,5)*elmas(5,j)/elmas(5,5)
     end do
  end do

  errnr = 0
  return

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

1000 return

end subroutine elem8
