      subroutine bpdcsto()

c     Unterprogramm berechnet b = B * p . 
c     Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
c
c     Copyright by Andreas Kemna 2009
c     
c     Andreas Kemna / Roland Martin                            10-Jun-2009
c
c     Letzte Aenderung   RM                                    30-Jul-2010
c
c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod
      USE modelmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      real            * 8     dum
      integer         * 4     i,j

c.....................................................................
      ALLOCATE (pvec2dc(manz),stat=errnr)
      IF (errnr /= 0) THEN
         fetxt = 'Error memory allocation pvec2dc in bpdcsto'
         errnr = 94
         RETURN
      END IF

c     A * p  berechnen (skaliert)
      do i=1,nanz
         apdc(i) = 0d0

         if (ldc) then
            do j=1,manz
               apdc(i) = apdc(i) + pvecdc(j)*sensdc(i,j)*fak(j)
            end do
         else if (lip) then
            do j=1,manz
               apdc(i) = apdc(i) + pvecdc(j)*dble(sens(i,j))*fak(j)
            end do
         end if
      end do

c     R^m * p  berechnen (skaliert)
caa   Abgeändert auf (4 Zeilen)
      do i=1,manz
         pvec2dc(i)=pvecdc(i)*fak(i)
      end do

      bvecdc = MATMUL(smatm,pvec2dc)

c     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
      do j=1,manz
         dum = 0d0

         if (ldc) then
            do i=1,nanz
               dum = dum + sensdc(i,j)*
     1              wmatd(i)*dble(wdfak(i))*apdc(i)
            end do
         else if (lip) then
            do i=1,nanz
               dum = dum + dble(sens(i,j))*
     1              wmatd(i)*dble(wdfak(i))*apdc(i)
            end do
         end if

         bvecdc(j) = dum + lam*bvecdc(j)
         bvecdc(j) = bvecdc(j)*fak(j)
      end do

      DEALLOCATE (pvec2dc)

      end
