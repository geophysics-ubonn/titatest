      subroutine broughtri()
c
c     Belegen der Leitfaehigkeit und zum Bestimmen der Rauhigkeit...
c     Fuer beliebige Triangulierung
c
c     Andreas Kemna                                            12-Apr-1996
c
c     Letzte Aenderung                                         29-Jul-2009
c
c.....................................................................
      USE alloci

      INCLUDE 'parmax.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'model.fin'
      INCLUDE 'sigma.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'

!.....................................................................

!     PROGRAMMINTERNE PARAMETER:

!     Hilfsvariablen
      integer         * 4     i,j,k,smaxs
      complex         * 16    cdum

!.....................................................................

      
!     Roughness bestimmen
      rough = 0d0
      smaxs=MAXVAL(selanz)

      DO i=1,manz
         cdum = dcmplx(0d0)
         DO j=1,nachbar(i,0)
            IF (nachbar(i,j)/=0)cdum=cdum+
     1           DCMPLX(smatm(i,j))*par(nachbar(i,j))
         END DO
         cdum = cdum + dcmplx(smatm(i,smaxs+1))*par(i)
         if (lip) then
            rough = rough + dimag(cdum)*dimag(par(i))
         else
            rough = rough + dble(cdum*dconjg(par(i)))
         end if
         
         
         
      END DO
!     Leitfaehigkeiten belegen
      do k=1,elanz
         j = mnr(k)
         sigma(k) = cdexp(par(j))
      end do
      

      return
      end
 
