      subroutine bptri(bvec,pvec)

c     Unterprogramm berechnet b = B * p .

c     Andreas Kemna                                            29-Feb-1996
c     Letzte Aenderung                                         29-Jul-2009
c...................................................................

      USE alloci

      INCLUDE 'parmax.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'model.fin'
      INCLUDE 'fem.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'elem.fin'

!.....................................................................

!     EIN-/AUSGABEPARAMETER:

!     Vektoren
      complex         * 16    bvec(mmax)
      complex         * 16    pvec(mmax)

!.....................................................................

!     PROGRAMMINTERNE PARAMETER:

!     Hilfsvektor
      complex         * 16    ap(nmax)

!     Hilfsvariablen
      complex         * 16    cdum
      integer         * 4     i,j,idum
      integer         * 4     smaxs ! maximale knotenanzahl

!.....................................................................
!     
!     A * p  berechnen (skaliert)

      smaxs=MAXVAL(selanz)

      do i=1,nanz
         ap(i) = dcmplx(0d0)

         do j=1,manz
            ap(i) = ap(i) + pvec(j)*sens(i,j)*dcmplx(fak(j))
         end do
      end do
      

!     R^m * p  berechnen (skaliert)
      DO i=1,manz
         cdum = dcmplx(0d0)
         DO j=1,nachbar(i,0)
            idum=nachbar(i,j)
            IF (idum/=0) cdum = cdum + pvec(idum)*
     1           DCMPLX(smatm(i,j))*fak(idum)
         END DO

         bvec(i) = cdum + pvec(i)*dcmplx(smatm(i,smaxs+1))*fak(i)

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
      

      return
      end

