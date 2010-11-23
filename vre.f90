subroutine vre()

!!!$     Fuehrt das Vorwaerts- und Rueckwaertseinsetzen mit der Cholesky-Links-
!!!$     dreiecksmatrix aus;
!!!$     'b' bleibt unveraendert, 'pot' ist Loesungsvektor.

!!!$     ( Vgl. Subroutine 'VRBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   24-Feb-1997

!!!$.....................................................................

  USE alloci
  USE femmod
  USE elemmod

  IMPLICIT none

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER (KIND=4)   :: idi,i0
  INTEGER (KIND=4)   :: m1,jlow
  COMPLEX(KIND(0D0)) ::   s

!!!$     Indexvariablen
  INTEGER (KIND=4)   ::     i,j

!!!$.....................................................................

  m1 = mb+1

  do i=1,sanz
     idi  = i*m1
     s    = b(i)
     i0   = idi-i
     jlow = max0(1,i-mb)

     do j=jlow,i-1
        s = s - a(i0+j)*pot(j)
     END do

     pot(i) = s / a(idi)
  END do

  do i=sanz,1,-1
     pot(i) = - pot(i) / a(idi)

     jlow = max0(1,i-mb)
     i0   = idi-i

     do j=jlow,i-1
        pot(j) = pot(j) + a(i0+j)*pot(i)
     END do

     idi = idi-m1
  END do

  return
end subroutine vre
