      SUBROUTINE batadcreg
c     
c     Unterprogramm berechnet ATC_d^-1A+lam*C_m
c     Fuer beliebige Triangulierung
c     
c     Andreas Kemna                                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   06-Nov-2009
c     
c.........................................................................
      USE alloci
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                      :: i,j,k,smaxs
!.....................................................................

c$$$  A^TC_d^-1A+lamC_m
      smaxs=MAXVAL(selanz)
      ALLOCATE (atadcreg(manz,manz),STAT=errnr)
      IF (errnr /= 0) THEN
         errnr = 97
         RETURN
      END IF

      IF (ltri == 0) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  atadcreg(i,j) = atadc(i,j) + lam * smatm(i,1)
                  IF (i+1 < manz) atadcreg(i+1,j) = atadc(i+1,j) + 
     1                 lam * smatm(i+1,2)
                  IF (i+nx < manz) atadcreg(i+1,j) = atadc(i+1,j) + 
     1                 lam * smatm(i+nx,3)
                  IF (i-1 > 1) atadcreg(i-1,j) = atadc(i-1,j) + 
     1                 lam * smatm(i-1,2)
                  IF (i-nx > 1) atadcreg(i-nx,j) = atadc(i-nx,j) + 
     1                 lam * smatm(i-nx,3)
               ELSE
                  atadcreg(i,j) = atadc(i,j) 
               END IF
            END DO
         END DO
      ELSE IF (ltri < 10 ) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  atadcreg(i,j) = atadc(i,j) + lam * smatm(i,smaxs+1)
                  DO k=1,nachbar(i,0)
                     IF (nachbar(i,k) /= 0) atadcreg(nachbar(i,k),j) = 
     1                    atadc(nachbar(i,k),j) + lam * smatm(i,k)
                  END DO
               ELSE
                  atadcreg(i,j) = atadc(i,j)
               END IF
            END DO
         END DO
      ELSE
         atadcreg = atadc + smatm    ! for full C_m..
      END IF

      PRINT*,MINVAL(atadcreg),MAXVAL(atadcreg)
      END
