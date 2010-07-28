      subroutine broughlma

c     Unterprogramm zum Belegen der Leitfaehigkeit und zum Bestimmen der
c     Rauhigkeit fuer Levenberg-Marquardt 
c
c     Copyright by Andreas Kemna 2010
c     
c     Andreas Kemna / Roland Martin                            24-Feb-2010
c
c     Letzte Aenderung   RM                                    24-Feb-2010
c
c.....................................................................

      USE alloci,only:smatm
      USE invmod
      USE sigmamod
      USE modelmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     i,j,k
      complex         * 16    cdum

c.....................................................................

c     Roughness bestimmen
      rough = 0d0

      DO j=1,manz
         if (.not.lprior) then

            cdum = DCMPLX(smatm(j,1)) * par(j)
            
            if (lip) then
               do i=1,manz
                  rough = rough + dimag(cdum)*dimag(par(i))
               end do
            else
               do i=1,manz
                  rough = rough + dble(cdum*dconjg(par(i)))
               end do
            end if
c     diff+<
         else
            
            cdum = DCMPLX(smatm(j,1)) * (par(j)-m0(j))

            if (lip) then
               do i=1,manz
                  rough = rough + dimag(cdum)*dimag(par(i)-m0(i))
               end do
            else
               do i=1,manz
                  rough = rough + dble(cdum*dconjg(par(i)-m0(i)))
               end do
            end if
         end if

      END DO

      end
