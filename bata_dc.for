      SUBROUTINE bata_dc(kanal)
c     
c     Unterprogramm berechnet ATC_d^-1A
c
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
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

      ALLOCATE(dig(manz))

      ata_dc = 0.
      DO k=1,manz
         write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1        'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
         DO j=k,manz ! fills upper triangle (k,j)
            DO i=1,nanz
               ata_dc(k,j) = ata_dc(k,j) + sensdc(i,k) * 
     1              sensdc(i,j) * wmatd(i) * DBLE(wdfak(i))
            END DO
            ata_dc(j,k) = ata_dc(k,j) ! fills lower triangle (k,j)
         END DO
         dig(k) = ata_dc(k,k)
      END DO

      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)

      errnr = 1
      OPEN (kanal,file=TRIM(fetxt)//'_re',
     1     status='replace',err=999)
      errnr = 4
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
