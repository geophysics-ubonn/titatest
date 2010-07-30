      subroutine bsmatmlma      ! levenberg-marquardt damping
c     
c     Unterprogramm belegt die Glaettungsmtrix
c     (hier Daempfungsmatrix (smatm))
c     
c     Copyright by Andreas Kemna                               2009
c     Erstellt von Roland Martin                               18-Dec-2009
c     Letzte Aenderung RM                                      Jul-2010
c     
c.........................................................................

      USE alloci
      USE femmod ! fuer ldc..
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
      INTEGER         :: i,l,k,j
!.....................................................................
      
      IF (.NOT.ALLOCATED (smatm))
     1     ALLOCATE (smatm(manz,1),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem smatm in bsmatmlma'
         errnr = 97
         RETURN
      END IF

      smatm = 0d0               ! initialize smatm

      IF (ltri==3) THEN
         
         PRINT*,'pure Levenberg (LA)'
         smatm = 1.0 ! Levenberg Damping

      ELSE

         PRINT*,' Levenberg-Marquardt (LMA)'
         IF (lip) THEN
            DO j=1,manz
               DO i=1,nanz
                  smatm(j,1) = smatm(j,1) + 
     1                 DBLE(sens(i,j))*DBLE(sens(i,j))* 
     1                 wmatd(i)*DBLE(wdfak(i))
               END DO
            END DO
         ELSE IF (ldc) THEN
            DO j=1,manz
               DO i=1,nanz
                  smatm(j,1) = smatm(j,1) + 
     1                 sensdc(i,j)*sensdc(i,j)* 
     1                 wmatd(i)*DBLE(wdfak(i))
               END DO
            END DO
         ELSE
            DO j=1,manz
               DO i=1,nanz
                  smatm(j,1) = smatm(j,1) + 
     1                 DCONJG(sens(i,j))*sens(i,j)* 
     1                 wmatd(i)*dble(wdfak(i)) ! wechselt automatisch zu
!     wmatdp bei lip
               END DO
            END DO
         ENDIF

      END IF

      END

