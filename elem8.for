        subroutine elem8(xk,yk,kelmas,kelve,kwert)

c Unterprogramm liefert die Elementmatrix 'kelmas(4,4)' und den Element-
c vektor 'kelve(4)' fuer ein zusammengesetztes Viereckelement
c (vier Teildreiecke mit linearem Ansatz) ( Elementtyp Nr.8 ).

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   26-Jan-1998
        
c.....................................................................

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'

c.....................................................................

c EIN-/AUSGABEPARAMETER:

c x-Koordinaten der Eckknotenpunkte
        real            * 8     xk(selmax)

c y-Koordinaten der Eckknotenpunkte
        real            * 8     yk(selmax)

c Elementmatrix nach der Kondensation
        real            * 8     kelmas(selmax,selmax)

c Elementvektor nach der Kondensation
        real            * 8     kelve(selmax)

c Wellenzahlwert
        real            * 8     kwert

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfskoordinaten
        real            * 8     xkdum(5),
     1                          ykdum(5)

c Beteiligte Knoten des zusammengesetzten Elements am Teildreieck
        integer         * 4     ik(3)

c Elementmatrizen vor der Kondensation
        real            * 8     elmas(5,5),
     1                          elmam(5,5)

c Elementvektor vor der Kondensation
        real            * 8     elve(5)

c Indexvariablen
        integer         * 4     i,j,k

c Grundelementmatrizen
        integer         * 4     s1(3,3),s2(3,3),s3(3,3),
     1                          s4(3,3)

c Grundelementvektor
        integer         * 4     sb(3)

c Hilfsvariablen
        real            * 8     x21,x31,y21,y31,
     1                          det,a,b,c,
     1                          a1,a2

c.....................................................................

        data s1/1,-1,0,-1,1,0,0,0,0/
        data s2/2,-1,-1,-1,0,1,-1,1,0/
        data s3/1,0,-1,0,0,0,-1,0,1/
        data s4/2,1,1,1,2,1,1,1,2/
        data sb/1,1,1/

c Elementmatrix bzw. -vektor des zusammengesetzten Elements gleich
c Null setzen
        do i=1,5
            elve(i) = 0d0

            do j=1,5
                elmas(i,j) = 0d0
                elmam(i,j) = 0d0
            end do
        end do

c Hilfskoordinaten bestimmen
        a1 = xk(3)-xk(1)
        a2 = xk(4)-xk(2)

c     Ggf. Fehlermeldung
        if (dabs(a1).le.1d-12.or.dabs(a2).le.1d-12) then
            fetxt = ' '
cak            write(fetxt(1:20),'(g20.5)') dabs(a1)
cak            write(fetxt(26:45),'(g20.5)') dabs(a2)
            errnr = 26
            goto 1000
        end if

        a1 = (yk(3)-yk(1))/a1
        a2 = (yk(4)-yk(2))/a2

c     Ggf. Fehlermeldung
        if (dabs(a1-a2).le.1d-12) then
            fetxt = ' '
cak            write(fetxt(1:20),'(g20.5)') dabs(a1-a2)
            errnr = 26
            goto 1000
        end if

        xkdum(5) = (a1*xk(1)-yk(1)-a2*xk(2)+yk(2))/(a1-a2)
        ykdum(5) = yk(1)+a1*(xkdum(5)-xk(1))

c Restlichen Koordinaten umspeichern
        do i=1,4
	      xkdum(i) = xk(i)
	      ykdum(i) = yk(i)
        end do

c Beitraege der Teilelemente aufaddieren
        do k=1,4

c         Beteiligte Knotenpunkte definieren
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

c         Ggf. Fehlermeldung
            if (det.le.1d-12) then
                fetxt = ' '
cak                write(fetxt(1:20),'(g20.5)') det
                errnr = 26
                goto 1000
            end if

            a =   (x31*x31 + y31*y31) / det
            b = - (x31*x21 + y31*y21) / det
            c =   (x21*x21 + y21*y21) / det

            do i=1,3
                elve(ik(i)) = elve(ik(i)) + det * dble(sb(i)) / 6d0

                do j=1,3
                    elmas(ik(i),ik(j)) = elmas(ik(i),ik(j)) +
     1                                   (a*dble(s1(i,j)) +
     1                                    b*dble(s2(i,j)) +
     1                                    c*dble(s3(i,j))) / 2d0
                    elmam(ik(i),ik(j)) = elmam(ik(i),ik(j)) +
     1                                   det * dble(s4(i,j)) / 2.4d1
                end do
            end do
        end do

c Kondensationsschritt
        do i=1,5
            do j=1,5
                elmas(i,j) = elmas(i,j) + kwert*kwert*elmam(i,j)
            end do
        end do

        do i=1,4
            kelve(i) = elve(i) - elmas(i,5)*elve(5)/elmas(5,5)

            do j=1,4
                kelmas(i,j) = elmas(i,j) -
     1                        elmas(i,5)*elmas(5,j)/elmas(5,5)
            end do
        end do

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

1000    return

        end
