      SUBROUTINE bata_reg(kanal)
c     
c     Unterprogramm berechnet ATC_d^-1A+lam*C_m
c     Fuer beliebige Triangulierung
c     
c     Copyright Andreas Kemna
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
c     
c.........................................................................
      USE alloci
      USE invmod
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
      INTEGER                       :: kanal ! io number
!     Hilfsvariablen 
      INTEGER                       :: i,j,k,smaxs
      REAL,DIMENSION(:),ALLOCATABLE :: dig
      REAL                          :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A+lamC_m
      
      errnr = 1
      OPEN(kanal,file=TRIM(fetxt),status='replace',err=999)
      errnr = 4

      ALLOCATE(dig(manz))

      smaxs=MAXVAL(selanz)

      IF (ltri == 0) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,j ! lower triangle
               
               IF (i == j) THEN
                  ata_reg(i,j) = ata(i,j) + DCMPLX(lam * smatm(i,1))
                  
                  IF (i+1 < manz) ata_reg(i+1,j) = ata(i+1,j) + 
     1                 DCMPLX(lam * smatm(i+1,2)) ! nebendiagonale in x richtung
                  IF (i+nx < manz) ata_reg(i+nx,j) = ata(i+nx,j) +
     1                 DCMPLX(lam * smatm(i+nx,3)) ! nebendiagonale in z richtung
               ELSE
                  ata_reg(i,j) = ata(i,j) ! only aTa
               END IF

               ata_reg(j,i) = ata_reg(i,j) ! upper triangle 

            END DO
            dig(j) = REAL(ata_reg(j,j))
         END DO

      ELSE IF (ltri == 3.OR.ltri == 4) THEN
         ata_reg = ata
         DO i=1,manz
            ata_reg(i,i) = ata(i,j) + DCMPLX(lam * smatm(i,1))
            dig(i) = REAL(ata_reg(i,i))
         END DO

      ELSE IF (ltri == 1.OR.ltri == 2.OR.
     1        (ltri > 4 .AND. ltri < 15)) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  ata_reg(i,j) = ata(i,j) + lam * smatm(i,smaxs+1)
                  DO k=1,nachbar(i,smaxs+1)
                     IF (nachbar(i,k) /= 0) ata_reg(nachbar(i,k),j) = 
     1                    ata(nachbar(i,k),j) + DCMPLX(lam * smatm(i,k))
                  END DO
               ELSE
                  ata_reg(i,j) = ata(i,j)
               END IF
            END DO
            dig(j) = REAL(ata_reg(j,j))
         END DO
      ELSE IF (ltri == 15) THEN
         ata_reg = ata + DCMPLX(lam * smatm) ! for full C_m..
         DO j=1,manz
            dig(j) = REAL(ata_reg(j,j))
         END DO
      END IF
      
c     write out real and imaginary part
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
