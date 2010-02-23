      SUBROUTINE bata_reg_dc(kanal)
c     
c     Unterprogramm berechnet ATC_d^-1A+lam*C_m
c     Fuer beliebige Triangulierung
c
c     Copyright Andreas Kemna
c     Erstellt von Roland Martin                               02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
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
      INTEGER                                      :: kanal ! io number
!     Hilfsvariablen 
      INTEGER                                      :: i,j,k,smaxs
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
      REAL(KIND(0D0))                              :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A+lamC_m
      
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ALLOCATE(dig(manz))
      smaxs=MAXVAL(selanz)

      IF (ltri == 0) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i > j) THEN
                  IF (i == j + 1) THEN
                     ata_reg_dc(i,j) = ata_dc(i,j) + 
     1                    lam * smatm(i,2) ! nebendiagonale in x richtung
                  ELSE IF (i == j + nx) THEN
                     ata_reg_dc(i,j) = ata_dc(i,j) +
     1                    lam * smatm(i,3) ! nebendiagonale in z richtung
                  ELSE
                     ata_reg_dc(i,j) = ata_dc(i,j)
                  END IF
               ELSE IF (i < j) THEN
                  ata_reg_dc(i,j) = ata_reg_dc(j,i) ! symmetry
               ELSE
                  ata_reg_dc(i,j) = ata_dc(i,j) + lam * smatm(i,1)
               END IF
            END DO
            dig(j) = ata_reg_dc(j,j)
         END DO
      ELSE IF (ltri < 10 ) THEN
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
            dig(j) = ata_reg_dc(j,j)
         END DO
      ELSE
         ata_reg_dc = ata_dc + smatm ! for full C_m..
         DO j=1,manz
            dig(j) = ata_reg_dc(j,j)
         END DO
      END IF

      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)

      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(dig(i)),LOG10(dig(i)/dig_max)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig)

      errnr = 0
 999  RETURN

      END
