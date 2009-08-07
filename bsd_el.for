      SUBROUTINE bsd_el
c     
c     Unterprogramm zum Bestimmen der kleinsten moeglichen 
c     Skalenlaenge zur stochastischen Regularisierung
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          29-Jul-2009
c     
c     Letzte Aenderung                                         07-Aug-2009
c     
c.....................................................................

      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin'        ! fuer nachbar, nrel etc. 
      INCLUDE 'model.fin'       ! fuer manz

c     PROGRAMMINTERNE PARAMETER:-------------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik,jk
c     Maximale knotenanzahl nicht entarteter Elemente
      INTEGER :: smaxs
c     Schwerpunktskoordinaten der Flaechenelemente
      REAL(KIND(0D0)) :: spx1,spx2,spy1,spy2
c     Abstand
      REAL(KIND(0D0)) :: r
c-----------------------------------------------------------------------

      smaxs = selanz(1)
      sd_el = 1.e10;
      
      DO i=1,elanz

         spx1=0.;spy1=0.
         DO ik=1,smaxs
            spx1 = spx1 + sx(snr(nrel(i,ik)))
            spy1 = spy1 + sy(snr(nrel(i,ik)))
         END DO
         spx1 = spx1/smaxs; spy1 = spy1/smaxs

         DO j=1,elanz
            IF (j==i) CYCLE

            spx2=0.;spy2=0.
            DO jk=1,smaxs
               spx2 = spx2 + sx(snr(nrel(j,jk)))
               spy2 = spy2 + sy(snr(nrel(j,jk)))
            END DO
            spx2 = spx2/smaxs; spy2 = spy2/smaxs

            r = SQRT((spx1-spx2)**2 + (spy1-spy2)**2)
            sd_el = MIN(sd_el,r)

         END DO                 ! inner loop j=1,elanz
      END DO                    ! outer loop i=1,elanz

      WRITE (*,'(A,F10.4)')'Minimalabstand:: ',sd_el

      END SUBROUTINE bsd_el
