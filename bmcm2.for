      SUBROUTINE bmcm2(kanal)
c     
c     Unterprogramm berechnet Modellkovarianz nach 
c     Fehlerfortpflanzung (Gubbins 2004)
c     MCM = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A * 
c           (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung und Complex
c
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                      02-Nov-2009
c     
c     Letzte Aenderung    RM                             30-Mar-2010
c     
c....................................................................

      USE alloci
      USE datmod
      USE invmod
      USE sigmamod
      USE modelmod

      IMPLICIT none

      INCLUDE 'err.fin'
!....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                     :: kanal
      INTEGER                                     :: i,j,k
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig
      COMPLEX(KIND(0D0))                          :: dum
      REAL(KIND(0D0))                             :: dig_min,dig_max,p
!....................................................................
!     vorher wurde schon 
!     ata_reg = (A^TC_d^-1A + C_m^-1)^-1 A^TC_d^-1A berechnet
!     cov_m = (A^TC_d^-1A + C_m^-1)^-1

      ata = MATMUL (ata_reg,cov_m)
      
      ALLOCATE (dig(manz))

      DO i=1,manz
         dig(i) = ata(i,i)
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
         WRITE (kanal,*)LOG10(SQRT(DBLE(dig(i)))),DBLE(dig(i))
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
         dum = dcmplx(1d0)/sigma(i)
         p = DBLE(1d3*datan2(dimag(dum),dble(dum)))
         WRITE (kanal,*)DIMAG(dig(i))*p,DIMAG(dig(i))
      END DO
      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
      CLOSE(kanal)


      DEALLOCATE (dig)

      errnr = 0
 999  RETURN
      
      END
