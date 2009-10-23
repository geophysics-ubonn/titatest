      subroutine brough()

c     Unterprogramm zum Belegen der Leitfaehigkeit und zum Bestimmen der
c     Rauhigkeit.

c     Andreas Kemna                                            12-Apr-1996
c     Letzte Aenderung   15-Jan-2001

c.....................................................................
      USE alloci 
      IMPLICIT none
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
c cdum describes (R^TR)m
c.....................................................................

c     Roughness bestimmen
      rough = 0d0

      do i=1,manz
         cdum = dcmplx(0d0)

c     diff+<
         if (.not.lprior) then
c     diff+>
            if (i.gt.1)
     1           cdum = dcmplx(smatm(i-1,2))*par(i-1)
            if (i.lt.manz)
     1           cdum = cdum + dcmplx(smatm(i,2))*par(i+1)
            if (i.gt.nx)
     1           cdum = cdum + dcmplx(smatm(i-nx,3))*par(i-nx)
            if (i.lt.manz-nx+1)
     1           cdum = cdum + dcmplx(smatm(i,3))*par(i+nx)

            cdum = cdum + dcmplx(smatm(i,1))*par(i)

            if (lip) then
               rough = rough + dimag(cdum)*dimag(par(i))
            else
               rough = rough + dble(cdum*dconjg(par(i)))
            end if
c     diff+<
         else
            if (i.gt.1) cdum = dcmplx(smatm(i-1,2)) *
     1           (par(i-1) - m0(i-1))
            if (i.lt.manz) cdum = cdum + dcmplx(smatm(i,2)) *
     1           (par(i+1) - m0(i+1))
            if (i.gt.nx) cdum = cdum + dcmplx(smatm(i-nx,3)) * 
     1           (par(i-nx) - m0(i-nx))
            if (i.lt.manz-nx+1) cdum = cdum + dcmplx(smatm(i,3)) * 
     1           (par(i+nx) - m0(i+nx))
            
            cdum = cdum + dcmplx(smatm(i,1))*(par(i)-m0(i))

            if (lip) then
               rough = rough + dimag(cdum) * dimag(par(i) - m0(i))
            else
               rough = rough + dble(cdum*dconjg(par(i) - m0(i)))
            end if
         end if
c     diff+>
      end do

c     Leitfaehigkeiten belegen
      do k=1,elanz
         j = mnr(k)
         sigma(k) = cdexp(par(j))
      end do

      return
      end
