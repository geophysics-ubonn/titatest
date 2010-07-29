      subroutine bpdctri()
c     
c     Unterprogramm berechnet b = B * p .
c     Fuer beliebige Triangulierung
c     
c     Copyright by Andreas Kemna         2009
c
c     Created by Roland Martin                            29-Jul-2009
c      
c     Last changes      RM                                   Jul-2010
c     
c.....................................................................

      USE alloci , ONLY : sens,sensdc,smatm,nachbar
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
      real            * 8     dum
      integer         * 4     i,j
!.....................................................................
!     WRITE(*,*) "BPDCtri",errnr 
!     A * p  berechnen (skaliert)

      do i=1,nanz
         apdc(i) = 0d0

         if (ldc) then
            do j=1,manz
               apdc(i) = apdc(i) + pvecdc(j) * sensdc(i,j) * fak(j)
            end do
         else if (lip) then
            do j=1,manz
               apdc(i) = apdc(i) + pvecdc(j) * dble(sens(i,j)) * fak(j)
            end do
         end if
      end do

!     R^m * p  berechnen (skaliert)
      DO i=1,manz
         dum = 0d0
         DO j=1,smaxs
            IF (nachbar(i,j) /= 0) dum = dum + pvecdc(nachbar(i,j))* 
     1           smatm(i,j) * fak(nachbar(i,j)) ! off diagonals
         END DO
!     main diagonal
         bvecdc(i) = dum + pvecdc(i) * smatm(i,smaxs+1) * fak(i) 
      END DO

!     A^h * R^d * A * p + l * R^m * p  berechnen (skaliert)
      do j=1,manz
         dum = 0d0

         if (ldc) then
            do i=1,nanz
               dum = dum + sensdc(i,j)* 
     1              wmatd(i) * dble(wdfak(i)) * apdc(i)
            end do
         else if (lip) then
            do i=1,nanz
               dum = dum + dble(sens(i,j))* 
     1              wmatd(i) * dble(wdfak(i)) * apdc(i)
            end do
         end if

         bvecdc(j) = dum + lam * bvecdc(j)

         bvecdc(j) = bvecdc(j) * fak(j)
      end do

      end
