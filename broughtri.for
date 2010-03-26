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
      IMPLICIT none
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
! cdum describes (R^TR)m
!.....................................................................
!     Roughness bestimmen
      rough = 0d0
      smaxs=MAXVAL(selanz)
      IF (.NOT.lprior) THEN
         DO i=1,manz
            cdum = dcmplx(0d0)
            DO j=1,smaxs
               IF (nachbar(i,j)/=0)cdum = cdum +
     1              DCMPLX(smatm(i,j))*par(nachbar(i,j))
            END DO
            cdum = cdum + dcmplx(smatm(i,smaxs+1))*par(i)
            if (lip) then
               rough = rough + dimag(cdum)*dimag(par(i))
            else
               rough = rough + dble(cdum*dconjg(par(i)))
            end if
         END DO
      ELSE 
         DO i=1,manz
            cdum = dcmplx(0d0)
            DO j=1,smaxs
               IF (nachbar(i,j)/=0)cdum=cdum + DCMPLX(smatm(i,j))*
     1              (par(nachbar(i,j)) - m0(nachbar(i,j)))
            END DO
            cdum = cdum + dcmplx(smatm(i,smaxs+1)) * 
     1           (par(i) - m0(i))
            if (lip) then
               rough = rough + dimag(cdum) * dimag(par(i) - m0(i))
            else
               rough = rough + dble(cdum * dconjg(par(i)) - m0(i))
            end if
         END DO
      END IF

      end
 
