      SUBROUTINE bmcm(kanal,ols)
c     
c     Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
c     MCM = (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung
c     
c     Copyright Andreas Kemna
c     Andreas Kemna / Roland Martin                      02-Nov-2009
c     
c     Letzte Aenderung    RM                             30-Mar-2010
c     
c....................................................................

      USE alloci
      USE tic_toc
      USE modelmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'err.fin'
!....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                     :: i,kanal,j,c1
      COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: dig,dig2
      REAL(KIND(0D0))                             :: dig_min,dig_max,p
      LOGICAL,INTENT(IN),OPTIONAL                 :: ols 
      CHARACTER(80)                               :: csz
!....................................................................

c$$$  invert A^TC_d^-1A + C_m^-1
      errnr = 1

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (dig(manz),dig2(manz)) ! help 

c     get time
      CALL TIC(c1)

      cov_m = ata_reg

      IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN
         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1        'Factorization...'
c$$$         CALL ZPOTRF('U',manz,cov_m,manz,errnr)
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
c$$$     1        'Inverting...'
c$$$         CALL ZPOTRI('U',manz,cov_m,manz,errnr)
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
c$$$     1        'Filling lower Cov...'
c$$$         DO i= 1,manz
c$$$            WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')
c$$$     1           ACHAR(13)//ACHAR(9)//ACHAR(9)//
c$$$     1           ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'
c$$$            DO j = i+1,manz
c$$$               cov_m(j,i)=cov_m(i,j)
c$$$            END DO
c$$$         END DO
         CALL CHOLZ(cov_m,dig,manz,errnr)
         IF (errnr /= 0) THEN
            PRINT*,'Zeile::',cov_m(abs(errnr),:)
            PRINT*,'Spalte::',cov_m(:,abs(errnr))
            errnr = 108
            RETURN
         END IF
         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1        'Inverting...'
         CALL LINVZ(cov_m,dig,manz)

         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1        'Filling lower Cov...'

         DO i= 1 , manz

            WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')
     1           ACHAR(13)//ACHAR(9)//ACHAR(9)//
     1           ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'

            DO j = 1 , i - 1

               cov_m(i,j) = cov_m(j,i)

            END DO
         END DO
      ELSE
         WRITE (*,'(a)',ADVANCE='no')
     1        'Inverting Matrix (Gauss elemination)'
         CALL gauss_cmplx(cov_m,manz,errnr)

         IF (errnr /= 0) THEN
            PRINT*,'Zeile::',cov_m(abs(errnr),:)
            PRINT*,'Spalte::',cov_m(:,abs(errnr))
            errnr = 108
            RETURN
         END IF

      END IF

      csz = 'solution time'
      CALL TOC(c1,csz)

      work = MATMUL(cov_m,ata_reg)

      DO i=1,manz
         dig(i) = cov_m(i,i)
         dig2(i) = work(i,i)
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
         WRITE (kanal,*)LOG10(SQRT(DBLE(dig(i)))),DBLE(dig2(i))
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
         p = DBLE(1d3*DATAN2(DIMAG(dig(i)),DBLE(dig(i))))
         WRITE (kanal,*)DIMAG(dig(i)),p
      END DO
      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min(Im):',dig_max,'/',dig_min
      CLOSE(kanal)


      DEALLOCATE (dig)
      IF (ALLOCATED(work)) DEALLOCATE (dig2,work)

      errnr = 0
 999  RETURN
      END
