      subroutine bpsto()

c     Unterprogramm berechnet b = B * p .
c     Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
c     Fuer komplexes Modell
c
c     Copyright by Andreas Kemna 2009
c     
c     Andreas Kemna / Roland Martin                            10-Jun-2009
c
c     Letzte Aenderung   RM                                    30-Jun-2009
c
c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod
      USE modelmod

      IMPLICIT none

      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'

c.....................................................................
c     PROGRAMMINTERNE PARAMETER:
c     Hilfsvariablen
      complex         * 16    cdum
      integer         * 4     i,j

c.....................................................................
c     A * p  berechnen (skaliert)
      do i=1,nanz
         ap(i) = dcmplx(0d0)
         do j=1,manz
            ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
         end do
      end do

c     R^m * p  berechnen (skaliert)
      do j=1,manz
         bvec(j) = 0.
         DO i = 1,manz
            bvec(j) = bvec(j) + pvec(i) * DCMPLX(smatm(i,j)) * 
     1           DCMPLX(fak(j))
         END DO
      end do

c     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
      do j=1,manz
         cdum = dcmplx(0d0)

         do i=1,nanz
            cdum = cdum + dconjg(sens(i,j))*
     1           dcmplx(wmatd(i)*dble(wdfak(i)))*ap(i)
         end do

         bvec(j) = cdum + dcmplx(lam)*bvec(j)
         bvec(j) = bvec(j)*dcmplx(fak(j))
      end do

      end
