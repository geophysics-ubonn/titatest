      SUBROUTINE bvariogram
c     
c     Unterprogramm zum Bestimmen von Variogrammen
c     
c     Copyright by Andreas Kemna
c     
c     Erste Version von Roland Martin                          19-Mar-2010
c     
c     Letzte Aenderung                                         19-Mar-2010
c     
c.....................................................................
      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin'        ! fuer nachbar, nrel etc. 
      INCLUDE 'inv.fin'         ! fuer par

c     PROGRAMMINTERNE PARAMETER:-------------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik,jk,ifp
c Variogramm counters
      INTEGER :: vcx,vix,vcy,viy,vcz,viz
c     Maximale knotenanzahl nicht entarteter Elemente
      INTEGER :: smaxs
c     Schwerpunktskoordinaten der Flaechenelemente
      REAL(KIND(0D0)) :: spx1,spx2,spy1,spy2
c     Variogramm abstaende
      REAL(KIND(0D0)) :: rij,hx,hy
c     Variogramm varianzen
      REAL(KIND(0D0)) :: varij
c     ESP Abstaende
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: variox,varioz,varior
c-----------------------------------------------------------------------

      ALLOCATE (variox(elanz,2),varioz(elanz,2),varior(elanz,2))

      smaxs = selanz(1)

      vcx = 0; vcy = 0; vcz = 0

      DO i=1,elanz

         spx1=0.;spy1=0.
         DO ik=1,smaxs
            spx1 = spx1 + sx(snr(nrel(i,ik)))
            spy1 = spy1 + sy(snr(nrel(i,ik)))
         END DO
         spx1 = spx1/smaxs; spy1 = spy1/smaxs

         DO j=1,elanz
            
            spx2=0.;spy2=0.
            DO jk=1,smaxs
               spx2 = spx2 + sx(snr(nrel(j,ik)))
               spy2 = spy2 + sy(snr(nrel(j,ik)))
            END DO
            spx2 = spx2/smaxs; spy2 = spy2/smaxs
            
            hx = (spx1-spx2)
            hy = (spy1-spy2)
            rij = SQRT(hx**2. + hy**2.)
            varij = (REAL(par(i)) - REAL(par(j)))**2.

            IF (vcx == 0) THEN
               variox(1,1) =  varij
               variox(1,2) = hx
            ELSE
               DO vix=1,vcx
c                  PRINT*,vix
                  IF ((hx - variox(vix,2))<EPSILON(hx)) THEN
                     variox(vix,1) = (variox(vix,1) * (vix - 1) + varij) 
     1                    * vix
                  ELSE IF (hx>variox(vix,2)) THEN
                     vcx = vcx + 1
                     variox(vcx,2) = hx
                     variox(vcx,1) = varij * vcx
                  ELSE
                     PRINT*,'check',i,j,vix,vcx
                  END IF

               END DO
            END IF

         END DO                 ! inner loop j=1,elanz


      END DO                    ! outer loop i=1,elanz

      PRINT*,'vcx:',vcx

      CALL get_unit(ifp)
      OPEN (ifp,FILE='inv.variox',STATUS='replace')
      WRITE (ifp,'(a)')'#    Distance    Variogram'
      DO i = 1,vcx
         WRITE (ifp,'(2(G10.3,2X))') variox(i,:)
      END DO
      CLOSE (ifp)
      DEALLOCATE (variox,varioz,varior)

      END 
