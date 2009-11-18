      SUBROUTINE bata_dc(kanal)
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
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: work
      REAL(KIND(0D0))                              :: min_work,max_work
!.....................................................................

c$$$  A^TC_d^-1A


      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4
      ALLOCATE(work(manz))
      ata_dc = 0.
      DO k=1,manz
         write(*,'(a,1X,F6.2,A)',advance='no')ACHAR(13)//
     1        'ATC_d^-1A/ ',REAL( k * (100./manz)),'%'
         DO j=1,manz
            DO i=1,nanz
               ata_dc(k,j) = ata_dc(k,j) + sensdc(i,k) * 
     1              sensdc(i,j) * wmatd(i) * DBLE(wdfak(i))
            END DO
            IF (k /= j) ata_dc(j,k) = ata_dc(k,j)
         END DO
         work(k) = ata_dc(k,k)
      END DO

      DO i=1,manz
         
      END DO

      min_work = MINVAL(work)
      max_work = MAXVAL(work)

      PRINT*,min_work,max_work
      work = work/max_work ! normierung
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(work(i)),work(i)*max_work
      END DO
      CLOSE(kanal)
      errnr = 0
 999  RETURN
      
      END
