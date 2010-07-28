      subroutine elem1()

c     Unterprogramm liefert die Elementmatrix 'elmam(2,2)' und den Element-
c     vektor 'elve(2)' fuer ein Randelement mit linearem Ansatz ( Element-
c     typ Nr.1 ).

c     Andreas Kemna                                            11-Oct-1993
c     Letzte Aenderung   22-Sep-1998

c.....................................................................
      
      USE elemmod,only:xk,yk,elmam,elve

      IMPLICIT none
      
c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Grundelementmatrix
      integer         * 4     s4(2,2)

c     Grundelementvektor
      integer         * 4     sb(2)

c     Hilfsvariablen
      real            * 8     x21,y21,l

c     Indexvariablen
      integer         * 4     i,j

c.....................................................................

      data sb/1,1/
      data s4/2,1,1,2/

      x21 = xk(2) - xk(1)
      y21 = yk(2) - yk(1)

      l = dsqrt(x21*x21 + y21*y21)

      do i=1,2
         elve(i) = l * dble(sb(i)) / 2d0

         do j=1,2
            elmam(i,j) = l * dble(s4(i,j)) / 6d0
         end do
      end do

      return
      end
