      subroutine bpdc(bvec,pvec)

c     Unterprogramm berechnet b = B * p .

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung   09-Jan-1998

c.....................................................................

      USE alloci
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Vektoren
      real            * 8     bvec(mmax)
      real            * 8     pvec(mmax)

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
      do i=1,manz
         dum = 0d0

         if (i.gt.1)
     1        dum = pvec(i-1)*smatm(i-1,2)*fak(i-1)
         if (i.lt.manz)
     1        dum = dum + pvec(i+1)*smatm(i,2)*fak(i+1)
         if (i.gt.nx)
     1        dum = dum + pvec(i-nx)*smatm(i-nx,3)*fak(i-nx)
         if (i.lt.manz-nx+1)
     1        dum = dum + pvec(i+nx)*smatm(i,3)*fak(i+nx)

         bvec(i) = dum + pvec(i)*smatm(i,1)*fak(i)
      end do

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
