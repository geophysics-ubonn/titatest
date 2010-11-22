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
      integer         * 4     idi,i0
      integer         * 4     m1,jlow
      complex         * 16    s

!!!$     Indexvariablen
      integer         * 4     i,j

!!!$.....................................................................

      m1 = mb+1

      do 20 i=1,sanz
         idi  = i*m1
         s    = b(i)
         i0   = idi-i
         jlow = max0(1,i-mb)

         do 10 j=jlow,i-1
            s = s - a(i0+j)*pot(j)
 10      continue

         pot(i) = s / a(idi)
 20   continue

      do 40 i=sanz,1,-1
         pot(i) = - pot(i) / a(idi)

         jlow = max0(1,i-mb)
         i0   = idi-i

         do 30 j=jlow,i-1
            pot(j) = pot(j) + a(i0+j)*pot(i)
 30      continue

         idi = idi-m1
 40   continue

      return
      end
