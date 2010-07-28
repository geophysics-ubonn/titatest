      SUBROUTINE bata(kanal)
c     
c     Unterprogramm berechnet ATA=A^T C_d^{-1} A (complex)
c     
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
c
c.........................................................................

      USE alloci
      USE datmod
      USE invmod
      USE modelmod
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                       :: kanal
      INTEGER                       :: i,j,k
      REAL,DIMENSION(:),ALLOCATABLE :: dig
      REAL                          :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A

      errnr = 1
      OPEN(kanal,file=TRIM(fetxt),status='replace',err=999)
      errnr = 4

      ALLOCATE(dig(manz))

      DO k=1,manz
         write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1        'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
         DO j=k,manz ! fills upper triangle (k,j) 
            DO i=1,nanz
               ata(k,j) = ata(k,j) + DCONJG(sens(i,k)) * 
     1              sens(i,j) * wmatd(i) * DBLE(wdfak(i))
            END DO
            ata(j,k) = ata(k,j) ! fills lower triangle (j,k) 
         END DO
         dig(k) = REAL(ata(k,k))
      END DO

c     write out 
      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(dig(i)/dig_max),dig(i)
      END DO
      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min
      CLOSE(kanal)

      DEALLOCATE (dig)

      errnr = 0
 999  RETURN

      END
