      SUBROUTINE bnachbar
      
c     Unterprogramm zum Bestimmen der Elementnachbarn
c     zur Realisierung der Triangulationsregularisierung
c
c     Andreas Kemna
ci.A. Roland Martin                                            29-Jul-2009
c     Letzte Aenderung   29-Jul-2009
c 
c.....................................................................
 

      IMPLICIT none

      INCLUDE 'parmax.fin' ! fuer die felddefinitionen in elem.fin
      INCLUDE 'elem.fin' ! fuer nachbar, nrel etc. 

c PROGRAMMINTERNE PARAMETER:

c Indexvariablen
      integer         * 4     i,j,k,i2,j2,k2
c Flächenelementzähler
      integer         * 4     ifl,ifl2
c Hilfsvariable
      integer         * 4     idum,idum2

      ifl=0
      nachbar=0
      DO i=1,typanz
         IF (typ(i)>10) EXIT    ! nur die flächenelemente die als erste kommen
         DO j=1,nelanz(i)       ! alle FL-elemente durchgehen
            ifl=ifl+1
            nachbar(ifl,0)=selanz(i) ! knotenpunkte speichern :(
            DO k=1,selanz(i)    ! alle knoten jeden Elements
               idum=nrel(ifl,k) ! zeigt auf den Knotenindex 
               ifl2=0
               DO j2=1,nelanz(i) ! suche nach gemeinsamen Knoten
                  ifl2=ifl2+1
                  DO k2=1,selanz(i)
                     idum2=nrel(ifl2,k2)
                     IF (idum==idum2) nachbar(ifl,k)=ifl2
                  END DO
               END DO
            END DO

            WRITE (*,*)'Nachbarn von element',ifl,'::',nachbar(ifl,:)

         END DO
         
      END DO

      END 
