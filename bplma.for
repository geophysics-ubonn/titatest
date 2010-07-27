      subroutine bplma()
c
c     Unterprogramm berechnet b = B * p . 
c     Angepasst an Levenberg-Marquardt-Daempfung
c
c     Copyright by Andreas Kemna 2010
c     
c     Andreas Kemna / Roland Martin                            24-Feb-2010
c
c     Letzte Aenderung   RM                                    24-Feb-2010
c
c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'

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

c     coaa R^m * p  berechnen (skaliert)

      do j=1,manz
         bvec(i)=pvec(i)*dcmplx(fak(i))*DCMPLX(smatm(i,1))
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
