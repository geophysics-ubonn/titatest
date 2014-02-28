subroutine elem3()

!!!$     Unterprogramm liefert die Elementmatrizen 'elmas(3,3)' und 'elmam(3,3)'
!!!$     sowie den Elementvektor 'elve(3)' fuer ein Dreieckelement mit linearem
!!!$     Ansatz ( Elementtyp Nr.3 ).

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   05-Nov-1997

!!!$.....................................................................
use alloci, only: prec
  USE elemmod,only:xk,yk,elmam,elmas,elve
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Grundelementmatrizen
  INTEGER (KIND = 4)  ::     s1(3,3),s2(3,3),s3(3,3),s4(3,3)

!!!$     Grundelementvektor
  INTEGER (KIND = 4)  ::     sb(3)

!!!$     Hilfsvariablen
  REAL (prec)    ::     x21,x31,y21,y31,det,a,b,c

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::     i,j

!!!$.....................................................................

  data s1/1,-1,0,-1,1,0,0,0,0/
  data s2/2,-1,-1,-1,0,1,-1,1,0/
  data s3/1,0,-1,0,0,0,-1,0,1/
  data s4/2,1,1,1,2,1,1,1,2/
  data sb/1,1,1/

  x21 = xk(2) - xk(1)
  x31 = xk(3) - xk(1)
  y21 = yk(2) - yk(1)
  y31 = yk(3) - yk(1)
  det = x21*y31 - x31*y21

!!!$     Ggf. Fehlermeldung
  if (det.le.0d0) then
     fetxt = TRIM(fetxt)//' hat evtl falsche Kontennummerierung'
     print*,'det,x21,y31,x31,y21',det,x21,y31,x31,y21
     errnr = 26
     return
  end if

  a =   (x31*x31 + y31*y31) / det
  b = - (x31*x21 + y31*y21) / det
  c =   (x21*x21 + y21*y21) / det

  do i=1,3
     elve(i) = det * REAL(sb(i)) / 6d0

     do j=1,3
        elmas(i,j) = (a*REAL(s1(i,j)) + &
             b*REAL(s2(i,j)) + c*REAL(s3(i,j))) / 2d0
        elmam(i,j) = det * REAL(s4(i,j)) / 2.4d1
     end do
  end do

  errnr = 0

  return
end subroutine elem3
