      subroutine bsmatmtri      !tri
c
c     Unterprogramm belegt die Rauhigkeitsmatrix....
c     fuer beliebige Triangulierung 
c     Angelehnt an R. Blaschek (2008)
c
c     Copyright by Andreas Kemna 2009
c     
c     Andreas Kemna / Roland Martin                            23-Jun-2009
c     
c     Letzte Aenderung   RM                                    Jul-2010
c
c.........................................................................
      USE alloci
      USE datmod
      USE invmod
      USE modelmod       ! mit nachbar und ldir
      USE elemmod
      USE errmod
      USE konvmod
      
      IMPLICIT none

!.....................................................................

!     PROGRAMMINTERNE PARAMETER:

!     Hilfsvariablen 
      REAL(KIND(0D0)) :: dum    ! dummy stores numbers
      INTEGER         :: i,l,k,ik,anz
      REAL(KIND(0D0)) :: edglen ! Kantenlaenge
      REAL(KIND(0D0)) :: dist   ! Abstand der Schwerpunkte
      REAL(KIND(0D0)) :: sp1(2),sp2(2) ! Schwerpunktkoordinaten
      REAL(KIND(0D0)) :: ang    !Winkel fuer anisotrope Glaettung
      REAL(KIND(0D0)) :: alfgeo !Anisotrope (geometrische) Glaettung
!.....................................................................
      
      IF (.NOT.ALLOCATED (smatm))
     1  ALLOCATE (smatm(manz,smaxs+1),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmtri'
         errnr = 97
         RETURN
      END IF

      smatm = 0D0               ! initialize smatm
      
      IF (elanz/=manz)PRINT*,'elanz/=manzSMATMTRI may be wrong'
      
      DO i=1,elanz
         
         sp1(1) = espx(i)
         sp1(2) = espy(i)
         
         DO k=1,smaxs           ! jedes flaechenele hat mind einen nachbarn
            
            ik = MOD(k,smaxs) + 1
            
            edglen = SQRT((sx(snr(nrel(i,k))) - 
     1           sx(snr(nrel(i,ik))))**2 +
     1           (sy(snr(nrel(i,k))) -
     1           sy(snr(nrel(i,ik))))**2) ! edge
            
            IF (nachbar(i,k)>0) THEN !nachbar existiert 
               
               sp2(1) = espx(nachbar(i,k)) ! schwerpunkt des nachbar elements
               sp2(2) = espy(nachbar(i,k))
!     Geometrischer Teil...

               dist = SQRT((sp1(1) - sp2(1))**2. + 
     1              (sp1(2) - sp2(2))**2.)
               
               ang = DATAN2((sp1(2) - sp2(2)),(sp1(1) - sp2(1))) !Winkel
               
               alfgeo = DSQRT((alfx*DCOS(ang))**2. + 
     1              (alfz*DSIN(ang))**2.)
               
               dum = edglen / dist * alfgeo
               
               smatm(i,k) = -dum
               
               smatm(i,smaxs+1) = smatm(i,smaxs+1) + dum !Hauptdiagonale

            END IF
            
         END DO
      END DO
      
      END

