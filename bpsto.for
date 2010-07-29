      subroutine bpsto()
c
c     Unterprogramm berechnet b = B * p .
c     Angepasst an die neue Regularisierungsmatrix 
c     (stoch. Kovarianzmatrix) fuer komplexes Modell
c
c     Copyright by Andreas Kemna 2009
c     
c     Created by Roland Martin                              10-Jun-2009
c
c     Last changes   RM                                     Jul-2010
c
c.....................................................................

      USE alloci , ONLY : sens,smatm
      USE femmod , ONLY : fak,ldc
      USE datmod , ONLY : nanz
      USE invmod , ONLY : lip,wmatd,wdfak
      USE cjgmod
      USE modelmod , ONLY : manz

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
      DO j=1,manz
         bvec(j) = 0.
         DO i = 1,manz
            bvec(j) = bvec(j) + pvec(i) * DCMPLX(smatm(i,j)) * 
     1           DCMPLX(fak(j))
         END DO
      END DO

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
