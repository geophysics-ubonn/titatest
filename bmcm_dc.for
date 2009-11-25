      SUBROUTINE bmcm_dc(kanal)
c     
c     Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
c     (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung
c     
c     Andreas Kemna                                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   06-Nov-2009
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
      INTEGER                                    :: i,kanal
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: ipiv
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: dig
      REAL(KIND(0D0))                            :: dig_min,dig_max
!....................................................................

c$$$  invert (A^TC_d^-1A + C_m^-1)

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

c$$$      ALLOCATE (work(manz,manz),STAT=errnr)
c$$$      IF (errnr/=0) THEN
c$$$         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
c$$$         errnr = 97
c$$$         RETURN
c$$$      END IF
c$$$      ALLOCATE (ipiv(manz),STAT=errnr)
c$$$      IF (errnr/=0) THEN
c$$$         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
c$$$         errnr = 97
c$$$         RETURN
c$$$      END IF
c$$$      
c$$$      work = ata_reg_dc         ! work is replaced by PLU decomposition
c$$$      
c$$$c$$$  building Right Hand Side (unit matrix)
c$$$      DO i=1,manz
c$$$         cov_m_dc(i,i) = 1.d0
c$$$      END DO
c$$$      
c$$$c$$$  Solving Linear System Ax=B -> B=A^-1
c$$$      WRITE (*,'(a)')'Solving Ax=B'
c$$$      CALL DGESV(manz,manz,work,manz,ipiv,cov_m_dc,manz,errnr)
c$$$      
c$$$      IF (errnr /= 0) THEN
c$$$         PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
c$$$         PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
c$$$         errnr = 108
c$$$         RETURN
c$$$      END IF
c$$$
c$$$      DEALLOCATE (work,ipiv)

      cov_m_dc = ata_reg_dc

      CALL gauss_dble(cov_m_dc,manz,errnr)
      
      IF (errnr /= 0) THEN
         fetxt = 'error matrix inverse not found'
         PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
         PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
         errnr = 108
         RETURN
      END IF
      
      ALLOCATE (dig(manz)) !prepare to write out main diagonal
      DO i=1,manz
         dig(i) = cov_m_dc(i,i)
      END DO
      
      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)
      
      PRINT*,dig_min,dig_max
      
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(SQRT(ABS(dig(i)))),dig(i)
      END DO
      CLOSE (kanal)
      DEALLOCATE (dig)

      errnr = 0
 999  RETURN

      END
