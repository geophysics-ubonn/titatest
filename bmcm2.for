      SUBROUTINE bmcm2(kanal)
c     
c     Unterprogramm berechnet Modellkovarianz nach Fehlerfortpflanzung
c     MCM = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung und Complex
c
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   23-Nov-2009
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
!     vorher wurde schon 
!     ata_reg = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!     cov_m = (A^TC_d^-1A + C_m^-1)^-1

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ata = MATMUL (ata_reg,cov_m)
      
      ALLOCATE (dig(manz))

      DO i=1,manz
         dig(i) = DBLE(ata(i,i))
      END DO
      
      dig_min = MINVAL(dig) 
      dig_max = MAXVAL(dig)
      
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(SQRT(ABS(dig(i)))),dig(i)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig)

      errnr = 0
 999  RETURN
      
      END
