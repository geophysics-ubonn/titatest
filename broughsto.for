      subroutine broughsto

c     Unterprogramm zum Belegen der Leitfaehigkeit und zum Bestimmen der
c     Rauhigkeit. 
c     Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).
c     
c     Copyright by Andreas Kemna 2009
c     
c     Andreas Kemna / Roland Martin                            10-Jun-2009
c     
c     Letzte Aenderung   RM                                    30-Jun-2009
c     
c.....................................................................

      USE sigmamod
      USE alloci,only:smatm
      USE invmod
      USE modelmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     i,j,k
      complex         * 16    cdum
      complex	      * 16    parh(manz)  

c     parh: Parameter-Hilfsvektor (R^TR)m bzw (C_m^-1)m

c.....................................................................

c     Roughness bestimmen

c     diff+<
      if (.not.lprior) then
c     diff+>

         parh = MATMUL(par(1:manz),DCMPLX(smatm))

         if (lip) then
            rough = DOT_PRODUCT(DIMAG(parh(1:manz)),
     1           DIMAG(par(1:manz)))
         else
            rough = DOT_PRODUCT(DBLE(parh(1:manz)),
     1           DCONJG(par(1:manz)))
         end if
         
c     diff+<
      else
         
         parh = MATMUL((par(1:manz)-m0(1:manz)),DCMPLX(smatm))
         IF (lip) THEN
            rough = DOT_PRODUCT(DIMAG(parh(1:manz)),
     1           DIMAG(par(1:manz)-m0(1:manz)))
         ELSE
            rough = DOT_PRODUCT(DBLE(parh(1:manz)),
     1           DCONJG(par(1:manz)-m0(1:manz)))
         END IF

      end if
c     diff+>
      end
