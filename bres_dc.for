      SUBROUTINE bres_dc(kanal)
c     
c     Unterprogramm berechnet Aufloesungsmatrix
c     RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A
c     Fuer beliebige Triangulierung
c     
c     Andreas Kemna                                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   23-Nov-2009
c     
c.........................................................................
      USE alloci
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                      :: i,kanal
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: dig
      REAL(KIND(0D0))                              :: dig_min,dig_max
!.....................................................................

c$$$  calc RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ata_reg_dc = MATMUL(cov_m_dc,ata_dc) ! that's it...

      ALLOCATE (dig(manz)) !prepare to write out main diagonal

      DO i=1,manz
         dig(i) = ata_reg_dc(i,i)
      END DO
      
      dig_min = MINVAL(ABS(dig))
      dig_max = MAXVAL(ABS(dig))
      
      PRINT*,dig_min,dig_max
      
      dig = dig/dig_max
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(ABS(dig(i))),dig(i)*dig_max
      END DO
      CLOSE (kanal)
      errnr = 0
      DEALLOCATE (dig)
 999  RETURN

      END
