      SUBROUTINE bata(kanal)
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
      INTEGER                                      :: kanal
      INTEGER                                      :: i,j,k
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
      REAL(KIND(0D0))                              :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A


      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4
      ALLOCATE(dig(manz))
      WRITE (kanal,*)manz
      DO k=1,manz
         write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1        'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
         DO j=1,manz
            DO i=1,nanz
               ata(k,j) = ata(k,j) + DCONJG(sens(i,k)) * 
     1              sens(i,j) * wmatd(i) * DBLE(wdfak(i))
            END DO
            IF (k /= j) ata(j,k) = ata(k,j)
         END DO
         dig(k) = DBLE(ata(k,k))
      END DO

      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)

      PRINT*,dig_min,dig_max
      dig = dig/dig_max ! normierung
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(dig(i)),dig(i)*dig_max
      END DO
      CLOSE(kanal)

      errnr = 0
 999  RETURN
      
      END
