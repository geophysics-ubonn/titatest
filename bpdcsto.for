      subroutine bpdcsto()
c
c     Unterprogramm berechnet b = B * p . 
c     Angepasst an die neue Regularisierungsmatrix
c     (stoch. Kovarianzmatrix)
c
c     Copyright by Andreas Kemna 2009
c     
c     Created by Roland Martin                             10-Jun-2009
c
c     Last changes        RM                                Jul-2010
c
c.....................................................................

      USE alloci , ONLY : sens,sensdc,smatm
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
      REAL(KIND(0D0)),ALLOCATABLE,DIMENSION(:) :: pvec2
      real            * 8     dum
      integer         * 4     i,j

c.....................................................................
c$$$      ALLOCATE (pvec2(manz),stat=errnr)
c$$$      IF (errnr /= 0) THEN
c$$$         fetxt = 'Error memory allocation pvec2 in bpdcsto'
c$$$         errnr = 94
c$$$         RETURN
c$$$      END IF

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
c$$$caa   Abgeändert auf (4 Zeilen)
c$$$      do i=1,manz
c$$$         pvec2(i)=pvecdc(i)*fak(i)
c$$$      end do
      do j = 1 , manz
         bvecdc(j) = 0.
         DO i = 1 , manz
            bvecdc(j) = bvecdc(j) + pvecdc(i) * smatm(i,j) * fak(i)
         END DO
      end do

c$$$
c$$$      bvecdc = MATMUL(smatm,pvec2)
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

      IF (ALLOCATED (pvec2)) DEALLOCATE (pvec2)

      end
