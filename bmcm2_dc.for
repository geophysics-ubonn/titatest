      SUBROUTINE bmcm2_dc(kanal)
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
!.....................................................................

c$$$  A^TC_d^-1A


      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ata_dc = MATMUL (ata_reg_dc,cov_m_dc)
      
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)log10(sqrt(abs(ata_dc(i,i)))),
     1        ata_dc(i,i)
      END DO

      CLOSE(kanal)
      errnr = 0

 999  RETURN
      
      END
