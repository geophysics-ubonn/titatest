      SUBROUTINE bmcm_dc(kanal,ols)
c     
c     Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
c     (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung
c     
c     Copyright Andreas Kemna 2009
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
c     
c.........................................................................

      USE alloci
      USE tic_toc
      USE modelmod
      USE errmod
      USE konvmod

      IMPLICIT none

!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                    :: i,kanal,j,c1
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dig,dig2
      REAL(KIND(0D0))                            :: dig_min,dig_max
      LOGICAL,INTENT(IN),OPTIONAL                :: ols
      CHARACTER(80)                              :: csz
!....................................................................

c$$$  invert (A^TC_d^-1A + C_m^-1)

      errnr = 1

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         fetxt = 'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (dig(manz),dig2(manz)) !prepare to write out main diagonal

c     get time
      CALL TIC(c1)

      cov_m_dc = ata_reg_dc
      
      IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN !default
         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1        'Factorization...'
c$$$         CALL DPOTRF('U',manz,cov_m_dc,manz,errnr)
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
c$$$     1           'Inverting...'
c$$$         CALL DPOTRI('U',manz,cov_m_dc,manz,errnr)
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
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
c$$$               cov_m_dc(j,i)=cov_m_dc(i,j)
c$$$            END DO
c$$$         END DO

         CALL CHOLD(cov_m_dc,dig,manz,errnr)
         IF (errnr /= 0) THEN
            fetxt='CHOLD mcm :: matrix not pos definite..'
            PRINT*,'Zeile(',abs(errnr),')'
            errnr = 108
            RETURN
         END IF
         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1           'Inverting...'

         CALL LINVD(cov_m_dc,dig,manz)

         WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//
     1        'Filling lower Cov...'

         DO i= 1 , manz

            WRITE (*,'(A,1X,F6.2,A)',ADVANCE='no')
     1           ACHAR(13)//ACHAR(9)//ACHAR(9)//
     1           ACHAR(9)//'/ ',REAL( i * (100./manz)),'%'

            DO j = 1 , i - 1

               cov_m_dc(i,j) = cov_m_dc(j,i)

            END DO
         END DO

      ELSE
         WRITE (*,'(a)',ADVANCE='no')
     1        'Inverting Matrix (Gauss elemination)'

         CALL gauss_dble(cov_m_dc,manz,errnr)
         
         IF (errnr /= 0) THEN
            fetxt = 'error matrix inverse not found'
            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
            errnr = 108
            RETURN
         END IF
      END IF
      
      csz = 'solution time'
      CALL TOC(c1,csz)

      work = MATMUL(cov_m_dc,ata_reg_dc)

      DO i=1,manz
         dig(i) = cov_m_dc(i,i)
         dig2(i) = work(i,i)
      END DO

      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)
      
      errnr = 1
      OPEN (kanal,file=TRIM(fetxt)//'_re',
     1     status='replace',err=999)
      errnr = 4

      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(SQRT(dig(i))),dig2(i)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig,dig2,work)

      errnr = 0
 999  RETURN

      END
