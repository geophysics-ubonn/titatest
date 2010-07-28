      SUBROUTINE bres(kanal)
c     
c     Unterprogramm berechnet Aufloesungsmatrix
c     Fuer beliebige Triangulierung und Complex
c     RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A
c     wobei (A^TC_d^-1A + C_m^-1) bereits invertiert wurde (cov_m)
c     
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                      02-Nov-2009
c     
c     Letzte Aenderung    RM                             30-Mar-2010
c     
c....................................................................

      USE alloci
      USE modelmod
      
      IMPLICIT none

      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'
!....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                     :: i,kanal
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig
      REAL(KIND(0D0))                             :: dig_min,dig_max
!     switch solv ordinary linear system or not^^
!....................................................................

c$$$  calculate  RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

      ata_reg = MATMUL(cov_m,ata)

      ALLOCATE (dig(manz))

      DO i=1,manz
         dig(i) = ata_reg(i,i)
      END DO
      
c     write out real and imaginary part
      errnr = 1
      OPEN(kanal,file=TRIM(fetxt)//'_re',
     1     status='replace',err=999)
      errnr = 4
      dig_min = MINVAL(DBLE(dig))
      dig_max = MAXVAL(DBLE(dig))
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(DBLE(dig(i))),DBLE(dig(i))
      END DO
      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min(Re):',dig_max,'/',dig_min
      CLOSE(kanal)

      errnr = 1
      OPEN(kanal,file=TRIM(fetxt)//'_im',
     1     status='replace',err=999)
      errnr = 4
      dig_min = MINVAL(DIMAG(dig))
      dig_max = MAXVAL(DIMAG(dig))
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(DIMAG(dig(i))),DIMAG(dig(i))
      END DO
      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
      CLOSE(kanal)

      DEALLOCATE (dig)

      errnr = 0

 999  RETURN

      END
