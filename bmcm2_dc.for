      SUBROUTINE bmcm2_dc(kanal)
c     
c     Unterprogramm berechnet Modellkovarianz nach Fehlerfortpflanzung
c     MCM = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung
c
c     Copyright Andreas Kemna/Roland Martin 2009
c     erstellt von Roland Martin                               02-Nov-2009
c     
c     Letzte Aenderung    RM                                   23-Nov-2009
c     
c.........................................................................
      USE alloci
      USE datmod
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'inv.fin'
      INCLUDE 'model.fin'
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
!     ata_reg_dc = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!     cov_m_dc = (A^TC_d^-1A + C_m^-1)^-1


      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ata_dc = MATMUL (ata_reg_dc,cov_m_dc)
      
      ALLOCATE (dig(manz))

      DO i=1,manz
         dig(i) = ata_dc(i,i)
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
