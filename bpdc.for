      subroutine bpdc()

c     Unterprogramm berechnet b = B * p .

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   09-Jan-1998

c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod
      USE cjgmod
      USE modelmod

      IMPLICIT none

      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      real            * 8     dum
      integer         * 4     i,j

c.....................................................................

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
      do i=1,manz
         dum = 0d0

         if (i.gt.1)
     1        dum = pvecdc(i-1)*smatm(i-1,2)*fak(i-1)
         if (i.lt.manz)
     1        dum = dum + pvecdc(i+1)*smatm(i,2)*fak(i+1)
         if (i.gt.nx)
     1        dum = dum + pvecdc(i-nx)*smatm(i-nx,3)*fak(i-nx)
         if (i.lt.manz-nx+1)
     1        dum = dum + pvecdc(i+nx)*smatm(i,3)*fak(i+nx)

         bvecdc(i) = dum + pvecdc(i)*smatm(i,1)*fak(i)
      end do

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

      end
