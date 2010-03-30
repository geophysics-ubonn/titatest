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
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE  :: dig
      REAL(KIND(0D0))                              :: dig_min,dig_max
!.....................................................................

c$$$  A^TC_d^-1A


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
         dig(k) = ata(k,k)
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
