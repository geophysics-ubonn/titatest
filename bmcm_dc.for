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

      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                    :: i,kanal
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: ipiv
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dig,dig2
      REAL(KIND(0D0))                            :: dig_min,dig_max
      LOGICAL,INTENT(IN),OPTIONAL                :: ols
!....................................................................

c$$$  invert (A^TC_d^-1A + C_m^-1)

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF

      CALL TIC

      IF (.NOT.PRESENT(ols).OR..NOT.ols) THEN !default
         ALLOCATE (ipiv(manz),STAT=errnr)
         IF (errnr/=0) THEN
            WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
            errnr = 97
            RETURN
         END IF
         
         work = ata_reg_dc
         cov_m_dc = 0.
c$$$  building Right Hand Side (unit matrix)
         DO i=1,manz
            cov_m_dc(i,i) = 1.d0
         END DO
         
c$$$  Solving Linear System Ax=B -> B=A^-1
         WRITE (*,'(a)',ADVANCE='no')'Solving Ax=B (DGESV)'
         CALL DGESV(manz,manz,work,manz,ipiv,cov_m_dc,manz,errnr)
         
         IF (errnr /= 0) THEN
            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
            errnr = 108
            RETURN
         END IF

      ELSE
         cov_m_dc = ata_reg_dc

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
c$$$         WRITE (*,'(a)',ADVANCE='no')
c$$$     1        'Factorization of LHS'
c$$$
c$$$         work = ata_reg_dc ! LU decomposition of LHS
c$$$
c$$$         work = 0.
c$$$         CALL chold(cov_m_dc,work,manz,errnr) ! -> L is stored in work
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$         WRITE (*,'(a)',ADVANCE='no')'inverting L'
c$$$
c$$$         cov_m_dc = 0.0
c$$$         CALL linv(work,cov_m_dc,manz)
c$$$  building Right Hand Side (unit matrix)
c$$$         cov_m_dc = 0.
c$$$
c$$$         DO i=1,manz
c$$$            cov_m_dc(i,i) = 1.d0
c$$$         END DO
c$$$
c$$$         CALL DPOTRF('U',manz,work,manz,errnr)
c$$$
c$$$         CALL DPOTRS('U',manz,manz,work,manz,cov_m_dc,
c$$$     1        manz,errnr)
c$$$
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',work(abs(errnr),:)
c$$$            PRINT*,'Spalte::',work(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$         
c$$$         
c$$$         WRITE (*,'(a)',ADVANCE='no')'inverting L'
c$$$         CALL MDPOTRI('U',manz,cov_m_dc,manz,errnr)
c$$$         IF (errnr /= 0) THEN
c$$$            PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
c$$$            PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
c$$$            errnr = 108
c$$$            RETURN
c$$$         END IF
c$$$
      END IF

      CALL TOC

      work = MATMUL(cov_m_dc,ata_reg_dc)

      ALLOCATE (dig(manz),dig2(manz)) !prepare to write out main diagonal
      DO i=1,manz
         dig(i) = cov_m_dc(i,i)
         dig2(i) = work(i,i)
      END DO
      
      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)
      
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(SQRT(ABS(dig(i)))),dig2(i)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig,dig2,work)
      IF (ALLOCATED(ipiv)) DEALLOCATE (ipiv)

      errnr = 0
 999  RETURN

      END
