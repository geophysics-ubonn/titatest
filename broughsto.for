      subroutine broughsto

c     Unterprogramm zum Belegen der Leitfaehigkeit und zum Bestimmen der
c     Rauhigkeit. Angepasst an die neue Regularisierungsmatrix (stoch. Kovarianzmatrix).

c     Andreas Kemna                                            12-Apr-1996
c     Letzte Aenderung   30-Jun-2009

c.....................................................................
      USE alloci,only:smatm

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'model.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Hilfsvariablen
      integer         * 4     i,j,k
      complex         * 16    cdum
      complex	      * 16    parh(manz)  

c     parh: Parameter-Hilfsvektor

c.....................................................................

c     Roughness bestimmen
      rough = 0d0

caa   do i=1,manz
caa   cdum = dcmplx(0d0)

c     diff+<
      if (.not.ldiff) then
c     diff+>


         parh=MATMUL(DCMPLX(smatm),par(1:manz))

         if (lip) then
            do i=1,manz
               rough = rough + dimag(parh(i))*dimag(par(i))
            end do
         else
            do i=1,manz
               rough = rough + dble(parh(i)*dconjg(par(i)))
            end do
         end if
c     diff+<
      else

         parh=MATMUL(dcmplx(smatm),(par(1:manz)-m0(i:manz)))

         if (lip) then
            do i=1,manz
               rough = rough + dimag(parh(i))*dimag(par(i)-m0(i))
            end do
         else
            do i=1,manz
               rough = rough + dble(parh(i)*dconjg(par(i)-m0(i)))
            end do
         end if
      end if
c     diff+>
caa   end do

c     Leitfaehigkeiten belegen
      do k=1,elanz
         j = mnr(k)
         sigma(k) = cdexp(par(j))
      end do

      return
      end
