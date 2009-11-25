      SUBROUTINE bres(kanal,ols)
c     
c     Unterprogramm berechnet Aufloesungsmatrix
c     Fuer beliebige Triangulierung
c     RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A
c     wobei (A^TC_d^-1A + C_m^-1) bereits invertiert wurde (cov_m)
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
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                         :: i,kanal
      COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE   :: work
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: ipiv
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE        :: dig
      REAL(KIND(0D0))                                 :: dig_min,dig_max
      LOGICAL,OPTIONAL                                :: ols 
!     switch solv ordinary linear system or not^^
!.....................................................................

c$$$  calculate  RES = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ata_reg = MATMUL(cov_m,ata)

      ALLOCATE (dig(manz))

      DO i=1,manz
         dig(i) = DBLE(ata_reg(i,i))
      END DO
      
      dig_min = MINVAL(dig) 
      dig_max = MAXVAL(dig)
      
      PRINT*,dig_min,dig_max
      
      dig = dig/dig_max
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(ABS(dig(i))),dig(i)*dig_max
      END DO
      CLOSE (kanal)
      DEALLOCATE (dig)
      errnr = 0
 999  RETURN

      END
