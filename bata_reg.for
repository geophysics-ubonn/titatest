      SUBROUTINE bata_reg(kanal)
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
      INTEGER                                      :: kanal ! io number
!     Hilfsvariablen 
      INTEGER                                      :: i,j,k,smaxs
!.....................................................................

c$$$  A^TC_d^-1A+lamC_m
      
      smaxs=MAXVAL(selanz)
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      WRITE (kanal,*)manz
      IF (ltri == 0) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i > j) THEN
                  IF (i == j + 1) THEN
                     ata_reg(i,j) = ata(i,j) + 
     1                    DCMPLX(lam * smatm(i,2)) 
!     nebendiagonale in x richtung
                  ELSE IF (i == j + nx) THEN
                     ata_reg(i,j) = ata(i,j) +
     1                    DCMPLX(lam * smatm(i,3)) 
!     nebendiagonale in z richtung
                  ELSE
                     ata_reg(i,j) = ata(i,j)
                  END IF
               ELSE IF (i < j) THEN
                  ata_reg(i,j) = ata_reg(j,i) ! symmetry
               ELSE
                  ata_reg(i,j) = ata(i,j) + DCMPLX(lam * smatm(i,1))
               END IF
c$$$  IF (i==j) THEN   ! sparse C_m
c$$$  ata_reg(i,j) = ata(i,j) + lam * smatm(i,1)
c$$$  IF (i+1 < manz) ata_reg(i+1,j) = ata(i+1,j) + 
c$$$  1                 lam * smatm(i+1,2)
c$$$  IF (i+nx < manz) ata_reg(i+nx,j) = ata(i+nx,j) + 
c$$$  1                 lam * smatm(i+nx,3)
c$$$  IF (i-1 > 1) ata_reg(i-1,j) = ata(i-1,j) + 
c$$$  1                 lam * smatm(i-1,2)
c$$$  IF (i-nx > 1) ata_reg(i-nx,j) = ata(i-nx,j) + 
c$$$  1                 lam * smatm(i-nx,3)
c$$$  ELSE
c$$$  ata_reg(i,j) = ata(i,j) 
c$$$  END IF
            END DO
            WRITE (kanal,*)ata_reg(j,j),j
         END DO
      ELSE IF (ltri < 10 ) THEN
         DO j=1,manz
            write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1           'ATC_d^-1A+reg/ ',REAL( j * (100./manz)),'%'
            DO i=1,manz
               IF (i==j) THEN   ! sparse C_m
                  ata_reg(i,j) = ata(i,j) + lam * smatm(i,smaxs+1)
                  DO k=1,nachbar(i,0)
                     IF (nachbar(i,k) /= 0) ata_reg(nachbar(i,k),j) = 
     1                    ata(nachbar(i,k),j) + DCMPLX(lam * smatm(i,k))
                  END DO
               ELSE
                  ata_reg(i,j) = ata(i,j)
               END IF
            END DO
            WRITE (kanal,*)ata_reg(j,j),j
         END DO
      ELSE
         ata_reg = ata + DCMPLX(lam * smatm) ! for full C_m..
         DO j=1,manz
            WRITE (kanal,*)ata_reg(j,j),j
         END DO
      END IF
      
      PRINT*,MINVAL(REAL(ata_reg)),MAXVAL(REAL(ata_reg))
      
      CLOSE (kanal)
      errnr = 0

 999  RETURN

      END
