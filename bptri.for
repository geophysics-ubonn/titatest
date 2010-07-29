      subroutine bptri()
c     
c     Unterprogramm berechnet b = B * p .
c     Fuer beliebige Triangulierung
c     
c     Copyright by Andreas Kemna         2009
c
c     Created by Roland Martin                           29-Jul-2009
c     
c     Last changes      RM                                  Jul-2010
c     
c...................................................................

      USE alloci , ONLY : sens,smatm,nachbar
      USE femmod , ONLY : fak,ldc
      USE datmod , ONLY : nanz
      USE invmod , ONLY : lip,wmatd,wdfak
      USE cjgmod
      USE modelmod , ONLY : manz
      USE elemmod, ONLY : smaxs

      IMPLICIT none

      INCLUDE 'konv.fin'

!.....................................................................

!     PROGRAMMINTERNE PARAMETER:

!     Hilfsvariablen
      complex         * 16    cdum
      integer         * 4     i,j,idum
!.....................................................................
!     
!     A * p  berechnen (skaliert)

      do i=1,nanz
         ap(i) = dcmplx(0d0)

         do j=1,manz
            ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
         end do
      end do
      

!     R^m * p  berechnen (skaliert)
      DO i=1,manz
         cdum = dcmplx(0d0)
         DO j=1,smaxs
            idum=nachbar(i,j)
            IF (idum/=0) cdum = cdum + pvec(idum)*
     1           DCMPLX(smatm(i,j)) * DCMPLX(fak(idum)) ! off diagonals
         END DO

         bvec(i) = cdum + pvec(i) * DCMPLX(smatm(i,smaxs+1)) * 
     1        DCMPLX(fak(i)) ! + main diagonal

      END DO
      

!     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
      do j=1,manz
         cdum = dcmplx(0d0)

         do i=1,nanz
            cdum = cdum + dconjg(sens(i,j))*dcmplx(wmatd(i)*
     1           dble(wdfak(i)))*ap(i)
         end do

         bvec(j) = cdum + dcmplx(lam)*bvec(j)
         
         bvec(j) = bvec(j)*dcmplx(fak(j))
      end do

      end

