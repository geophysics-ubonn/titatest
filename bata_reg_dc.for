      SUBROUTINE bata_reg_dc(kanal)
c     
c     Unterprogramm berechnet ATC_d^-1A+lam*C_m
c     Fuer beliebige Triangulierung
c     
c     Copyright Andreas Kemna
c     Erstellt von Roland Martin                           02-Nov-2009
c     
c     Letzte Aenderung    RM                               31-Mar-2010
c     
c.........................................................................

      USE alloci
      USE invmod
      USE modelmod
      USE elemmod
      
      IMPLICIT none

      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
      INTEGER                        :: kanal ! io number
!     Hilfsvariablen 
      INTEGER                        :: i,j,k
      REAL,DIMENSION(:),ALLOCATABLE  :: dig
      REAL                           :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A+lamC_m
      
      errnr = 1
      open(kanal,file=TRIM(fetxt),status='replace',err=999)
      errnr = 4

      ALLOCATE(dig(manz))
      
      IF (ltri == 0) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,j            ! lower triangle
               
               IF (i == j) THEN
                  ata_reg_dc(i,j) = ata_dc(i,j) + lam * smatm(i,1)
                  
                  IF (i+1 < manz) ata_reg_dc(i+1,j) = ata_dc(i+1,j) + 
     1                 lam * smatm(i+1,2) ! nebendiagonale in x richtung
                  IF (i+nx < manz) ata_reg_dc(i+nx,j) = ata_dc(i+nx,j) +
     1                 lam * smatm(i+nx,3) ! nebendiagonale in z richtung
               ELSE
                  ata_reg_dc(i,j) = ata_dc(i,j) ! only aTa
               END IF

               ata_reg_dc(j,i) = ata_reg_dc(i,j) ! upper triangle 

            END DO
            
            dig(j) = REAL(ata_reg_dc(j,j))
         END DO

      ELSE IF (ltri == 3.OR.ltri == 4) THEN
         ata_reg_dc = ata_dc
         DO i=1,manz
            dig (i) = ata_dc(i,i) + lam * smatm(i,1)
            ata_reg_dc(i,i) = dig(i)
         END DO

      ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1        (ltri > 4 .AND. ltri < 15)) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  ata_reg_dc(i,j) = ata_dc(i,j) + lam * smatm(i,smaxs+1)
                  DO k=1,nachbar(i,smaxs+1)
                     IF (nachbar(i,k) /= 0) ata_reg_dc(nachbar(i,k),j) = 
     1                    ata_dc(nachbar(i,k),j) + lam * smatm(i,k)
                  END DO
               ELSE
                  ata_reg_dc(i,j) = ata_dc(i,j)
               END IF
            END DO
            dig(j) = REAL(ata_reg_dc(j,j))
         END DO

      ELSE IF (ltri == 15) THEN
         ata_reg_dc = ata_dc + smatm ! for full C_m..w
         DO j=1,manz
            dig(j) = REAL(ata_reg_dc(j,j))
         END DO
      END IF

      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)

      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(dig(i)),dig(i)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig)

      errnr = 0
 999  RETURN

      END
