      SUBROUTINE batadc
c     
c     Unterprogramm berechnet ATC_d^-1A
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
      INCLUDE 'model.fin'
      INCLUDE 'dat.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                      :: i,j,k
!.....................................................................

c$$$  A^TC_d^-1A

      ALLOCATE (atadc(manz,manz),STAT=errnr)
      IF (errnr /= 0) THEN
         errnr = 97
         RETURN
      END IF
      atadc = 0D0
      DO k=1,manz
         write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1        'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
         DO j=1,manz
            DO i=1,nanz
               atadc(k,j) = atadc(k,j) + sensdc(i,k) * 
     1              sensdc(i,j) * wmatd(i) * DBLE(wdfak(i))
            END DO
            IF (k /= j) atadc(j,k) = atadc(k,j)
         END DO
      END DO

      PRINT*,MINVAL(atadc),MAXVAL(atadc)

      END
