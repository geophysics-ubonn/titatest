      SUBROUTINE bmcmdc
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
      INTEGER                                      :: i
      REAL(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE   :: work
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE     :: ipiv
!.....................................................................

c$$$  solve (A^TC_d^-1A + C_m^-1) x = B

      ALLOCATE (atadcreg_inv(manz,manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem MCM_1 in bmcmdc'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (work(manz,manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz))
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF

c$$$  building Right Hand Side (unit matrix)
      DO i=1,manz
         atadcreg_inv(i,i) = 1.0
      END DO
      work = atadcreg ! work is replaced by PLU decomposition
c$$$  Solving Linear System Ax=B -> B=A^-1
      WRITE (*,'(a)')ACHAR(9)//'Solving Ax=B'
      CALL DGESV(manz,manz,work,manz,ipiv,atadcreg_inv,manz,errnr)
      
      IF (errnr /= 0) THEN
         PRINT*,'Zeile::',atadcreg_inv(abs(errnr),:)
         PRINT*,'Spalte::',atadcreg_inv(:,abs(errnr))
         errnr = 108
         RETURN
      END IF
      
      DEALLOCATE (work,ipiv)

      END
