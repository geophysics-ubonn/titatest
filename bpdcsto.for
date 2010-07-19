      subroutine bpdcsto(bvec,pvec)

c     Unterprogramm berechnet b = B * p . 
c     Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
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

      IMPLICIT none
      INCLUDE 'parmax.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Vektoren
      real            * 8     bvec(mmax)
      real            * 8     pvec(mmax)
      real            * 8     pvecdum(manz)

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvektor
      real            * 8     ap(nmax)

c     Hilfsvariablen
      real            * 8     dum
      integer         * 4     i,j

c.....................................................................

c     A * p  berechnen (skaliert)
      do i=1,nanz
         ap(i) = 0d0

         if (ldc) then
            do j=1,manz
               ap(i) = ap(i) + pvec(j)*sensdc(i,j)*fak(j)
            end do
         else if (lip) then
            do j=1,manz
               ap(i) = ap(i) + pvec(j)*dble(sens(i,j))*fak(j)
            end do
         end if
      end do

c     R^m * p  berechnen (skaliert)
caa   Abgeändert auf (4 Zeilen)
      do i=1,manz
         pvecdum(i)=pvec(i)*fak(i)
c     PRINT*,'pvec bp::',i,pvec(i),fak(i)
      end do

      bvec(1:manz)=MATMUL(smatm,pvecdum(1:manz))

c     PRINT*,'bvec cg::',bvec(1:manz)


c     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
      do j=1,manz
         dum = 0d0

         if (ldc) then
            do i=1,nanz
               dum = dum + sensdc(i,j)*
     1              wmatd(i)*dble(wdfak(i))*ap(i)
            end do
         else if (lip) then
            do i=1,nanz
               dum = dum + dble(sens(i,j))*
     1              wmatd(i)*dble(wdfak(i))*ap(i)
            end do
         end if

         bvec(j) = dum + lam*bvec(j)
         bvec(j) = bvec(j)*fak(j)
      end do

      return
      end
