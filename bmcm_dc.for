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
      INTEGER                                      :: i,kanal
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE   :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: ipiv
!.....................................................................

c$$$  solve (A^TC_d^-1A + C_m^-1) x = B

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF


      WRITE (kanal,*)manz

      work = ata_reg_dc ! work is replaced by PLU decomposition

c$$$  building Right Hand Side (unit matrix)
      DO i=1,manz
         cov_m_dc(i,i) = 1.0
      END DO

c$$$  Solving Linear System Ax=B -> B=A^-1
      WRITE (*,'(a)')ACHAR(9)//'Solving Ax=B'
      CALL DGESV(manz,manz,work,manz,ipiv,cov_m_dc,manz,errnr)
      
      IF (errnr /= 0) THEN
         PRINT*,'Zeile::',cov_m_dc(abs(errnr),:)
         PRINT*,'Spalte::',cov_m_dc(:,abs(errnr))
         errnr = 108
         RETURN
      END IF
      
      DEALLOCATE (work,ipiv)
      
      DO i=1,manz
         WRITE (kanal,*)log10(sqrt(abs(cov_m_dc(i,i)))),
     1        cov_m_dc(i,i)
      END DO
      
      CLOSE (kanal)
      errnr = 0
 999  RETURN

      END
