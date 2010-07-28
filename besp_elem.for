      SUBROUTINE besp_elem
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

      USE alloci
      USE elemmod        ! fuer nachbar, nrel etc. 

      IMPLICIT none

      INCLUDE 'parmax.fin'      ! fuer die felddefinitionen in elem.fin

c     PROGRAMMINTERNE PARAMETER:-------------------------------------------
c     Indexvariablen
      INTEGER :: i,j,ik,jk,fp
c     Maximale knotenanzahl nicht entarteter Elemente
      INTEGER :: smaxs
c     Schwerpunktskoordinaten der Flaechenelemente
      REAL(KIND(0D0)) :: spx1,spx2,spy1,spy2,ax,ay,ar
c     ESP Abstaende
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: abst
c-----------------------------------------------------------------------

      IF (.NOT.ALLOCATED(abst)) ALLOCATE (abst(elanz))

      smaxs = selanz(1)

      grid_min = 10**5.; grid_max = 0.
      grid_minx = 10**5.; grid_maxx = 0.
      grid_miny = 10**5.; grid_maxy = 0.

      DO i=1,elanz

         spx1=0.;spy1=0.
         DO ik=1,smaxs
            spx1 = spx1 + sx(snr(nrel(i,ik)))
            spy1 = spy1 + sy(snr(nrel(i,ik)))
         END DO
         spx1 = spx1/smaxs; spy1 = spy1/smaxs

         DO ik=1,smaxs
            
            IF (nachbar(i,ik)==0) CYCLE
            
            spx2=0.;spy2=0.
            DO jk=1,smaxs
               spx2 = spx2 + sx(snr(nrel(nachbar(i,ik),jk)))
               spy2 = spy2 + sy(snr(nrel(nachbar(i,ik),jk)))
            END DO
            spx2 = spx2/smaxs; spy2 = spy2/smaxs
            
            abst(i) = SQRT((spx1-spx2)**2 + (spy1-spy2)**2)
            
         END DO                 ! inner loop ik=1,smaxs

         DO ik=1,elanz
            
            IF (i==ik) CYCLE
            
            spx2=0.;spy2=0.
            DO jk=1,smaxs
               spx2 = spx2 + sx(snr(nrel(ik,jk)))
               spy2 = spy2 + sy(snr(nrel(ik,jk)))
            END DO
            spx2 = spx2/smaxs; spy2 = spy2/smaxs
            
            ax = ABS(spx1-spx2)
            ay = ABS(spy1-spy2)
            ar = SQRT(ax**2. + ay**2.)

            grid_min = MIN(grid_min,ar)
            grid_max = MAX(grid_max,ar)
            grid_minx = MAX(MIN(grid_minx,ax),grid_min)
            grid_maxx = MAX(grid_maxx,ax)
            grid_miny = MAX(MIN(grid_miny,ay),grid_min)
            grid_maxy = MAX(grid_maxy,ay)

         END DO

      END DO                    ! outer loop i=1,elanz
      
c     maximaler wert aus der Menge der Nachbarmittelpunkte
      esp_max = MAXVAL(abst) 
      esp_min = MINVAL(abst)

      esp_mit = SUM(abst)/elanz

      esp_std = 0.
      DO i=1,elanz
         esp_std = esp_std + SQRT((abst(i) - esp_mit)**2.)
      END DO
      esp_std = esp_std / MAX(1, elanz - 1)
      
      CALL MDIAN1(abst,elanz,esp_med)

      CALL get_unit(fp)
      OPEN (fp,FILE='inv.gstat',STATUS='replace')
      WRITE (fp,'(a/)')'Grid statistics:'
      WRITE (fp,'(20X,A,I10)')'Gridcells:'//ACHAR(9),elanz
      WRITE (fp,'(20X,A,2F10.4)')'ESP Min/Max:'//ACHAR(9),
     1     esp_min,esp_max
      WRITE (fp,'(20X,A,2F10.4)')'GRID Min/Max:'//ACHAR(9),
     1     grid_min,grid_max
      WRITE (fp,'(20X,A,2F10.4)')'GRID-x Min/Max:'//ACHAR(9),
     1     grid_minx,grid_maxx
      WRITE (fp,'(20X,A,2F10.4)')'GRID-y Min/Max:'//ACHAR(9),
     1     grid_miny,grid_maxy
      WRITE (fp,'(20X,A,3F10.4)')'Mean/Median/Var:'//ACHAR(9),
     1     esp_mit,esp_med,esp_std
      CLOSE (fp)

      IF (ALLOCATED(abst)) DEALLOCATE (abst)

      END SUBROUTINE besp_elem
      
