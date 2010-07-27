      subroutine blam0()

c     Unterprogramm zum Bestimmen des Start-Regularisierungsparameters.

c     Andreas Kemna                                            20-Feb-1997
c     Letzte Aenderung   07-Mar-2003
      
c.....................................................................

      USE alloci
      USE femmod
      USE datmod
      USE invmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'model.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      complex         * 16    cdum
      real            * 8     dum

c     Indexvariablen
      integer         * 4     i,j,k

c.....................................................................

c     Start-Regularisierungsparameter bestimmen
      IF (nz<0) THEN
         IF (nz<-1) lammax = -REAL(nz)
         IF (nz==-1) lammax = MAX(REAL(manz),REAL(nanz))
         WRITE (*,'(2x,a,F12.1)')'taking easy lam_0 ',lammax
         RETURN
      END IF

      lammax = 0d0

      if (ldc) then
         do j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'blam0/ ',REAL( j * (100./manz)),'%'
            dum = 0d0

            do i=1,nanz
               do k=1,manz
                  dum = dum + sensdc(i,j)*sensdc(i,k)*
     1                 wmatd(i)*dble(wdfak(i))
               end do
            end do

            lammax = lammax + dabs(dum)
         end do

      else if (lip) then

         do j=1,manz
            dum = 0d0

            do i=1,nanz
               do k=1,manz
                  dum = dum + dble(sens(i,j))*dble(sens(i,k))*
     1                 wmatd(i)*dble(wdfak(i))
               end do
            end do

            lammax = lammax + dabs(dum)
         end do

      else

         do j=1,manz
            cdum = dcmplx(0d0)

            do i=1,nanz
               do k=1,manz
                  cdum = cdum + dconjg(sens(i,j))*sens(i,k)*
     1                 dcmplx(wmatd(i)*dble(wdfak(i)))
               end do
            end do

            lammax = lammax + cdabs(cdum)
         end do

      end if

      lammax = lammax/dble(manz)

      lammax = lammax * 2d0/(alfx+alfz)
c     ak Default
      lammax = lammax * 5d0

c     ak Synthetic Example (JoH)
c     ak        lammax = lammax * 1d1

c     ak MinFrac
c     ak        lammax = lammax * 5d1

c     ak Test
c     ak        lammax = lammax * 1d1

c     ak AAC
c     ak        lammax = lammax * 5d0

      return
      end
